"""Portable prediction certificate, without optimization or project imports."""
from pathlib import Path
from fractions import Fraction as F
import argparse,copy,importlib.util,json
HERE=Path(__file__).resolve().parent;STABILITY=HERE.parent/'illumination-testing-stability';CORE=HERE.parent/'illumination-testing'
spec=importlib.util.spec_from_file_location('stability',STABILITY/'verify_stability.py');s=importlib.util.module_from_spec(spec);spec.loader.exec_module(s)
ps=importlib.util.spec_from_file_location('prediction',HERE/'predict.py');predictor=importlib.util.module_from_spec(ps);ps.loader.exec_module(predictor)
require=s.require

def prepare(p):
    require(s.v.digest(CORE/'proof.json')==p['coreProofSHA256'],'Core binding')
    require(s.v.digest(STABILITY/'proof.json')==p['stabilityProofSHA256'],'Stability proof binding')
    require(s.v.digest(STABILITY/'verify_stability.py')==p['stabilityVerifierSHA256'],'Stability verifier binding')
    sp=s.v.read(STABILITY/'proof.json');data=s.prepare(sp);s.check_claims(sp,data)
    core=s.v.read(CORE/'proof.json');poly=s.v.kernels(CORE)
    white=[sum((a+b/2+c/3 for a,b,c in row),F(0)) for row in poly]
    raw={c['id']:c for c in core['duals']}
    chosen=sorted(data['certs'],key=lambda c:c['coordinate'])
    V=[list(map(F,raw[c['id']]['v'])) for c in chosen]
    return dict(data=data,sp=sp,V=V,white=white,certs=chosen)

def check(p,context):
    require(p['schema']=='illumination-affine-prediction-v1','Schema')
    require(p['assumptions']==dict(actualPassivity=True,whiteCalibrationProvided=True,
        whiteCalibrationErrorCharged=True,independentBoxMeasurementNoise=True,
        physicalSpectrumOutputRequired=False,sharedSPDFixedPredictorClaimed=False),'Calibration and information assumptions')
    d=context['data'];sp=context['sp'];certs=context['certs'];V=context['V'];white=context['white']
    require(list(map(F,p['nominalTestWhite']))==white[:15] and list(map(F,p['nominalTargetWhite']))==white[15:],'Nominal white integrals')
    require([[F(t) for t in row] for row in p['V']]==V,'Coefficient identity')
    epsilon=F(p['epsilon']);require(2*epsilon==d['delta'],'Noise versus pair tolerance')
    require(F(p['predictionLimit'])==d['tau']/2,'Prediction target')
    mu,rho=F(sp['muAxis']),F(sp['rhoAxis'])
    for c in certs:
        for total_error,distance in [(0,0),(mu,0),(0,rho)]:
            require(c['base']/2+(1+c['N'])*total_error/2+c['U']*distance<=d['tau']/2,'Calibration region')
    for w in d['witnesses']:
        require(mu*w['amp']/2<=epsilon and (w['target']-mu*w['amp'])/2>d['tau']/2,'Common observation lower bound')
    examples=[]
    for e in p['examples']:
        m,ga,r=map(F,(e['mu'],e['gamma'],e['rho']))
        require(min(m,ga,r)>=0,'Negative error allowance')
        upper=max(c['base']/2+(1+c['N'])*(m+ga)/2+c['U']*r for c in certs)
        lower=min((w['target']-m*w['amp'])/2 for w in d['witnesses'])
        require(upper<=F(e['upper'])<=d['tau']/2 and lower>=F(e['lower'])>d['tau']/2,'Reported prediction bounds')
        examples.append(dict(upper=str(upper),lower=str(lower)))
    grey=dict(y=[str(w/2) for w in white[:15]],testWhite=list(map(str,white[:15])),targetWhite=list(map(str,white[15:])))
    require(predictor.predict(p,grey)==[w/2 for w in white[15:]],'Exact grey evaluation')
    return dict(status='PASS',continuousDuals=9,exactWitnesses=5,excludedSubsets=31,
        wholeCalibrationRegionChecked=True,examples=examples,exactGreyPrediction=True,optimizerCalls=0)

def faults(p,context):
    bad=[]
    q=copy.deepcopy(p);q['V'][0][0]=str(F(q['V'][0][0])+F(1,100));bad.append(('altered coefficient',q))
    q=copy.deepcopy(p);q['assumptions']['whiteCalibrationProvided']=False;bad.append(('missing white side information',q))
    q=copy.deepcopy(p);q['assumptions']['whiteCalibrationErrorCharged']=False;bad.append(('uncharged calibration error',q))
    q=copy.deepcopy(p);q['epsilon']=str(2*F(q['epsilon']));bad.append(('doubled measurement noise',q))
    q=copy.deepcopy(p);q['examples'][0]['upper']='0';bad.append(('understated prediction error',q))
    records=[]
    for name,q in bad:
        try:check(q,context)
        except ValueError as exc:records.append(dict(case=name,rejected=True,reason=str(exc)))
        else:raise ValueError('Invalid certificate accepted: '+name)
    return records

def main():
    a=argparse.ArgumentParser();a.add_argument('--output');a.add_argument('--check-faults',action='store_true');args=a.parse_args()
    p=s.v.read(HERE/'proof.json');context=prepare(p);out=check(p,context)
    if args.check_faults:out['faultCases']=faults(p,context)
    text=json.dumps(out,indent=2)+'\n'
    if args.output:Path(args.output).write_text(text)
    else:print(text)
if __name__=='__main__':main()
