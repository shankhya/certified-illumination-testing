"""Verify the complete continuous frontier plus nine exact dual blends."""
from pathlib import Path
from fractions import Fraction as F
import argparse,copy,importlib.util,json,time
HERE=Path(__file__).resolve().parent;CORE=HERE.parent/'illumination-testing'
spec=importlib.util.spec_from_file_location('core',CORE/'verify.py')
v=importlib.util.module_from_spec(spec);spec.loader.exec_module(v)
require=v.require

def reconstructed(recipe,raw):
    a,b=[raw[name] for name in recipe['parents']];theta=F(recipe['theta'])
    require(theta in [F(1,4),F(1,2),F(3,4)],'Undeclared blend weight')
    require(a['coordinate']==b['coordinate']==recipe['coordinate'],'Mixed target coordinates')
    require(F(a['h'])==F(b['h'])==1,'Mixed flux meshes')
    record={k:recipe[k] for k in ['id','coordinate','R','P','N','a','b']}
    record.update(h='1',left='0',right='1/500',
        v=[str((1-theta)*F(x)+theta*F(y)) for x,y in zip(a['v'],b['v'])],
        flux=[str((1-theta)*F(x)+theta*F(y)) for x,y in zip(a['flux'],b['flux'])])
    return record

def updated_geometry(p,original,duals,lines):
    combined=dict(original,safePieces=p['safePieces'],maximumGap=p['maximumGap'],maximumGapDelta=p['maximumGapDelta'])
    result=v.verify_envelopes(combined,duals,lines)
    require(F(p['maximumGap'])<=F(p['reportedGapUpper']),'Reported gap rounding')
    require(p['uniformGapTargetMet']==result['uniformGapTargetMet'],'False precision success flag')
    return result

def faults(p,original,raw,poly,duals,lines):
    checks=[]
    def reject(name,operation,expected):
        try:operation()
        except v.InvalidCertificate as error:
            require(expected in str(error),'Unexpected rejection cause: '+str(error))
            checks.append(dict(name=name,status='REJECTED_AS_EXPECTED',reason=str(error)))
        else:raise RuntimeError('Invalid refinement accepted: '+name)
    bad=copy.deepcopy(p['recipes'][0]);bad['theta']='3/4'
    reject('altered blend weight with old bound',lambda:v.verify_dual(reconstructed(bad,raw),poly),'Dual test norm')
    bad=copy.deepcopy(p['recipes'][0]);bad['parents'][1]='R160-exact-h100.json:c7'
    reject('mixed target coordinates',lambda:reconstructed(bad,raw),'Mixed target coordinates')
    badp=copy.deepcopy(p);badp['maximumGap']='0'
    reject('understated whole-interval gap',lambda:updated_geometry(badp,original,duals,lines),'Maximum continuum gap identity')
    badp=copy.deepcopy(p);badp['uniformGapTargetMet']=True
    reject('false target-met report',lambda:updated_geometry(badp,original,duals,lines),'False precision success flag')
    return checks

def verify(include_faults=False):
    start=time.monotonic();p=v.read(HERE/'proof.json');original=v.read(CORE/'proof.json')
    require(p['schema']=='illumination-frontier-refinement-v1','Unknown refinement schema')
    require(v.digest(CORE/'proof.json')==p['coreProofSHA256'],'Core proof changed')
    require(v.digest(CORE/'verify.py')==p['coreVerifierSHA256'],'Core verifier changed')
    v.check_sources(CORE,original);poly=v.kernels(CORE);duals=[];raw={}
    for index,c in enumerate(original['duals'],1):
        duals.append(v.verify_dual(c,poly));raw[c['id']]=c
        if index%9==0:print('Original continuous duals checked: '+str(index),flush=True)
    base_count=len(duals);require(base_count==81 and len(p['recipes'])==9,'Frozen certificate counts')
    for recipe in p['recipes']:
        require(recipe['id'] not in raw,'Duplicate refined certificate')
        c=reconstructed(recipe,raw);duals.append(v.verify_dual(c,poly));raw[c['id']]=c
    weights={F(h):v.nodal_integrals(poly,F(h)) for h in {w['h'] for w in original['witnesses']}}
    witnesses={w['id']:v.verify_witness(w,weights[F(w['h'])]) for w in original['witnesses']}
    require(len(witnesses)==32,'Frozen witness count')
    lines,points=v.verify_failure_lines(original,witnesses)
    v.verify_envelopes(original,duals[:base_count],lines)
    frontier=updated_geometry(p,original,duals,lines);count=v.verify_count(original,duals,witnesses)
    out=dict(status='PASS',optimizerCalls=0,dualCount=len(duals),newBlendedDuals=9,
        exactWitnessCount=len(witnesses),strictFrontierPoints=len(points),frontier=frontier,
        testCount=count,coreUnchanged=True,elapsedSeconds=round(time.monotonic()-start,3))
    if include_faults:out['faultCases']=faults(p,original,raw,poly,duals,lines)
    return out

if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--output',type=Path);parser.add_argument('--check-faults',action='store_true');args=parser.parse_args()
    result=verify(args.check_faults)
    if args.output:args.output.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2))
