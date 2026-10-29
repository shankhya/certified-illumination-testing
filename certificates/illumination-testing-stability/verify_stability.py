"""Portable exact verification of residual-aware and shared-SPD stability."""
from pathlib import Path
from fractions import Fraction as F
from itertools import combinations
import argparse, importlib.util, json
HERE=Path(__file__).resolve().parent;CORE=HERE.parent/'illumination-testing'
spec=importlib.util.spec_from_file_location('core_verifier',CORE/'verify.py')
v=importlib.util.module_from_spec(spec);spec.loader.exec_module(v)
require=v.require

def check_configuration(p):
    require(p['schema']=='illumination-stability-v1','Unknown stability schema')
    require(v.digest(CORE/'proof.json')==p['coreProofSHA256'],'Core proof binding')
    require(v.digest(CORE/'verify.py')==p['coreVerifierSHA256'],'Core verifier binding')
    require(p['assumptions']==dict(actualPassivity=True,observerFixed=True,
        sharedSPDWithinXYZ=True,normalization='actual-white-Y',
        continuousRelativeShapeBound=True,unrestrictedGainPositive=True), 'Stability assumptions changed')
    require(tuple(map(F,(p['g'],p['delta'],p['tau'])))==(F(71,2000),F(1,10000),F(1,50)), 'Operating condition')

def abs_quadratics(row):
    low=high=F(0)
    for a,b,c in row:
        lo,hi=v.absolute_cubic((a,a+b/3,a+2*b/3+c/3,a+b+c),F(1),F(1,10**13)/len(row))
        low+=lo;high+=hi
    require(high-low<=F(1,10**12),'Quadratic norm enclosure too wide')
    return low,high

def absolute_witness_integral(row,nodes):
    require(len(nodes)==len(row)+1,'Witness integral dimensions')
    total=F(0)
    for j,(a,b,c) in enumerate(row):
        d=nodes[j];z=nodes[j+1]-d;cuts=[F(0),F(1)]
        if z and 0<-d/z<1:cuts.insert(1,-d/z)
        power=(a*d,b*d+a*z,c*d+b*z,c*z)
        for lo,hi in zip(cuts,cuts[1:]):
            total+=abs(sum((co*(hi**(k+1)-lo**(k+1))/(k+1) for k,co in enumerate(power)),F(0)))
    return total

def prepare(p):
    check_configuration(p);original=v.read(CORE/'proof.json');v.check_sources(CORE,original)
    poly=v.kernels(CORE);white=[sum((a+b/2+c/3 for a,b,c in row),F(0)) for row in poly]
    raw={c['id']:c for c in original['duals']};certs=[]
    g,delta,tau=map(F,(p['g'],p['delta'],p['tau']))
    for chosen in p['certificates']:
        source=raw[chosen['id']];checked=v.verify_dual(source,poly)
        j=checked['coordinate'];weights=list(map(F,source['v']))
        u=[tuple(poly[14+j][cell][power]-sum((weights[i]*poly[i][cell][power] for i in range(15)),F(0)) for power in range(3)) for cell in range(400)]
        _,U=abs_quadratics(u);S=white[14+j]
        for light in range(5):
            block=[tuple(sum((weights[3*light+c]*poly[3*light+c][cell][power] for c in range(3)),F(0)) for power in range(3)) for cell in range(400)]
            _,norm=abs_quadratics(block);S+=norm
        N=checked['N'];base=checked['R']+2*g*checked['P']+delta*N
        certs.append(dict(id=chosen['id'],coordinate=j,U=U,S=S,N=N,base=base,
            margin=tau-base,A=1+N,C=tau+delta*N+S))
    require(sorted(c['coordinate'] for c in certs)==list(range(1,10)),'Target coordinate coverage')
    raw_w={w['id']:w for w in original['witnesses']};nodal=v.nodal_integrals(poly,F(1));witnesses=[]
    for name in p['witnessIds']:
        source=raw_w[name];checked=v.verify_witness(source,nodal)
        require(F(source['h'])==1 and checked['g']<=g and checked['delta']==0,'Witness class and exact matching')
        nodes=list(map(F,source['nodes']));amp=max(map(abs,nodes))
        j=max(range(9),key=lambda k:abs(checked['targets'][k]));target=abs(checked['targets'][j])
        rows=[3*(i-1)+c for i in checked['tests'] for c in range(3)]
        I=[absolute_witness_integral(poly[i],nodes) for i in rows]
        Iq=absolute_witness_integral(poly[15+j],nodes)
        witnesses.append(dict(id=name,tests=checked['tests'],amp=amp,target=target,I=I,Iq=Iq))
    require(len(witnesses)==5,'Five witness pairs required')
    return dict(certs=certs,witnesses=witnesses,g=g,delta=delta,tau=tau)

def norm_bound(data,mu,rho):
    return max(c['base']+c['A']*mu+2*c['U']*rho for c in data['certs'])

def shape_bound(data,kappa,rho):
    require(0<=kappa<1,'Invalid white denominator')
    return max((c['base']+kappa*(data['delta']*c['N']+c['S'])+2*rho*c['U'])/(1-kappa) for c in data['certs'])

def certify_norm_point(data,mu,rho):
    require(mu>=0 and rho>=0 and norm_bound(data,mu,rho)<=data['tau'],'Point outside norm guarantee')

def check_claims(p,data):
    check_configuration(p)
    certs=data['certs'];witnesses=data['witnesses'];delta=data['delta'];tau=data['tau']
    require([c['id'] for c in certs]==[c['id'] for c in p['certificates']],'Certificate selection changed')
    for c,declared in zip(certs,p['certificates']):
        require(0<c['U']<=F(declared['Uupper']),'Residual norm understated')
        require(0<c['S']<=F(declared['Supper']),'Shared block norm understated')
    mu,rho,kappa=map(F,(p['muAxis'],p['rhoAxis'],p['kappaAxis']))
    require(mu>0 and rho>0 and 0<kappa<1,'Nonpositive triangle')
    require(norm_bound(data,mu,F(0))<=tau and norm_bound(data,F(0),rho)<=tau,'Norm triangle unsafe')
    require(shape_bound(data,kappa,F(0))<=tau and shape_bound(data,F(0),rho)<=tau,'Shape triangle unsafe')
    # An affine inequality holds on the entire triangle once checked at its vertices.
    for w in witnesses:
        require(mu*w['amp']<=delta and w['target']-mu*w['amp']>tau,'Norm witness persistence')
        require(max(kappa*I/(1-kappa) for I in w['I'])<=delta,'Shape witness acceptance')
        require((w['target']-kappa*w['Iq'])/(1+kappa)>tau,'Shape witness strict failure')
    count=0
    for n in range(5):
        for T in combinations(range(1,6),n):
            require(any(set(T)<=set(w['tests']) for w in witnesses),'Smaller protocol uncovered');count+=1
    require(count==31,'Protocol coverage count')
    require(F(p['reportedMuAxis'])<=mu and F(p['reportedRhoAxis'])<=rho and F(p['reportedKappaAxis'])<=kappa,'Axis rounding')
    checked_examples=[]
    for item in p['examples']:
        er=F(item['rho']);e=F(item['error'])
        if item['kind']=='kernel':
            upper=norm_bound(data,e,er);lower=min(w['target']-w['amp']*e for w in witnesses)
            tests=max(w['amp']*e for w in witnesses)
        else:
            require(item['kind']=='shape','Unknown example uncertainty type')
            upper=shape_bound(data,e,er);lower=min((w['target']-e*w['Iq'])/(1+e) for w in witnesses)
            tests=max(e*I/(1-e) for w in witnesses for I in w['I'])
        require(upper<=F(item['upper'])<=tau and lower>=F(item['lower'])>tau and tests<=delta,'Example report or persistence')
        checked_examples.append(dict(kind=item['kind'],error=float(e),rho=float(er),upper=float(upper),lower=float(lower)))
    return dict(status='PASS',optimizerCalls=0,continuousDualsRechecked=9,exactWitnessesReintegrated=5,
        excludedProtocols=count,minimumConstituentTests=5,muAxis=float(mu),rhoAxis=float(rho),
        kappaAxis=float(kappa),normAndShapeTrianglesVerified=True,examples=checked_examples)

def controls():
    lo,hi=abs_quadratics([(-F(1,2),F(1),F(0))]);require(lo==hi==F(1,4),'Residual analytic integral')
    require(absolute_witness_integral([(F(1),F(0),F(0))],[F(1),-F(1)])==F(1,2),'P1 absolute integral')
    require(absolute_witness_integral([(F(0),F(0),F(1))],[F(1),-F(1)])==F(3,16),'Cubic sign split')
    k=F(1,10);actual=(1-(1-k))/(1+(1-k));require(actual==k/(2-k)>k/2,'Denominator analytic control')
    _,norm=abs_quadratics([(F(0),F(0),F(0))]);require(norm==0,'Shared cancelling block')
    return 5

def verify():
    p=v.read(HERE/'proof.json');data=prepare(p);out=check_claims(p,data)
    out['analyticControls']=controls();return out

if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--output',type=Path);args=parser.parse_args()
    result=verify()
    if args.output:args.output.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2))
