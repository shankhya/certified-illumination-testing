"""Analytic controls and meaningful invalid-certificate rejection checks."""
from pathlib import Path
from fractions import Fraction as F
import argparse, copy, importlib.util, json

ROOT=Path(__file__).resolve().parent
spec=importlib.util.spec_from_file_location('certificate_verifier',ROOT/'verify.py')
v=importlib.util.module_from_spec(spec);spec.loader.exec_module(v)

def main():
    p=v.read(ROOT/'proof.json');poly=v.kernels(ROOT);out=[]
    def rejected(name,operation,expected):
        try:operation()
        except v.InvalidCertificate as error:
            v.require(expected in str(error),'Wrong rejection cause: '+str(error))
            out.append(dict(name=name,status='REJECTED_AS_EXPECTED',reason=str(error)))
        else:raise RuntimeError('Invalid certificate accepted: '+name)
    for beta,truth in [([F(2)]*4,F(2)),([0,F(1,3),F(1,3),0],F(1,6)),([0,-F(1,3),-F(1,3),0],F(1,6))]:
        lo,hi=v.absolute_cubic(beta,F(1),F(1,10**12));v.require(lo==hi==truth,'Analytic fixed-sign integral')
    primitive=lambda t:t**4/4-t**3/3+3*t*t/32
    truth=2*primitive(F(1,4))-2*primitive(F(3,4))+primitive(F(1))
    lo,hi=v.absolute_cubic([0,F(1,16),-F(5,24),F(3,16)],F(1),F(1,10**12))
    v.require(lo<=truth<=hi and hi-lo<=F(1,10**12),'Analytic sign-changing cubic enclosure')
    bad=copy.deepcopy(p);first=next(iter(bad['sourceData']));bad['sourceData'][first]['sha256']='0'*64
    rejected('source hash corruption',lambda:v.check_sources(ROOT,bad),'Source hash mismatch')
    dual=copy.deepcopy(p['duals'][0]);dual['flux'][0]='1/1000000'
    rejected('nonzero boundary flux',lambda:v.verify_dual(dual,poly),'endpoint must be zero')
    dual2=copy.deepcopy(p['duals'][0]);dual2['R']=str(F(dual2['R'])+F(1,10**12))
    rejected('changed residual cost',lambda:v.verify_dual(dual2,poly),'residual identity failed')
    original=next(w for w in p['witnesses'] if w['id'].startswith('R162:'))
    weights=v.nodal_integrals(poly,F(original['h']))
    valid=v.verify_witness(original,weights)
    badw=copy.deepcopy(original);badw['nodes'][0]='2'
    rejected('nonpassive witness',lambda:v.verify_witness(badw,weights),'passive amplitude')
    badw2=copy.deepcopy(original);d=list(map(F,badw2['nodes']));k=min(range(1,len(d)-1),key=lambda k:abs(d[k]));d[k]+=F(1,10**12)
    badw2['nodes']=list(map(str,d));badw2['g']=str(max(abs(y-x)/2 for x,y in zip(d,d[1:])))
    rejected('broken exact metamer',lambda:v.verify_witness(badw2,weights),'Exact four-light matching failed')
    # Supply trusted points only for this isolated interpolation-geometry check.
    mock={}
    for point in p['points']:
        w=next(w for w in p['witnesses'] if w['id']==point['witnessId'])
        mock[w['id']]=dict(tests=w['tests'],targets=list(map(F,w['targets'])),g=F(w['g']),delta=max(map(abs,map(F,w['testResponses']))))
    mixed=copy.deepcopy(p);line=next(line for line in mixed['failureLines'] if len(line['pointIds'])==2)
    other=next(point for point in p['points'] if point['coordinate']!=line['coordinate'])
    line['pointIds'][1]=other['id']
    rejected('mixed-coordinate witness interpolation',lambda:v.verify_failure_lines(mixed,mock),'Mixed-coordinate interpolation')
    return dict(status='PASS',optimizerCalls=0,analyticCases=4,signChangingIntegralExact=str(truth),faultCases=out,validExactMetamerControl=True)

if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--output',type=Path);args=parser.parse_args()
    result=main()
    if args.output:args.output.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2))
