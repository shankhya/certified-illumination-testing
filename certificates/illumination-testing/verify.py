"""Exact, optimizer-free verification of the bundled continuous certificates.

Python 3.10+ standard library only. No project imports or saved PASS flags.
The accompanying README states the mathematics and limits of this check.
"""
from fractions import Fraction as F
from pathlib import Path
from itertools import combinations
from datetime import datetime, timezone
import argparse, csv, hashlib, json, time

TAU, END = F(1,50), F(1,500)

class InvalidCertificate(ValueError):
    pass

def require(condition, message):
    if not condition:
        raise InvalidCertificate(message)

def read(path):
    return json.loads(path.read_text(encoding='utf-8-sig'))

def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()

def check_sources(root, proof):
    require(proof['schema']=='illumination-testing-certificate-v1', 'Unknown schema')
    require(proof['model']['wavelengthNm']==[380,780], 'Wavelength domain changed')
    require(proof['model']['tests']==['A','D65','FL2','FL7','FL11'], 'Test library changed')
    require(proof['model']['targets']==['FL3','FL8','FL12'], 'Deployment library changed')
    require(F(proof['model']['tau'])==TAU, 'Mismatch target changed')
    require(list(map(F,proof['model']['toleranceInterval']))==[0,END], 'Tolerance domain changed')
    for name, meta in proof['sourceData'].items():
        require(digest(root/'data'/name)==meta['sha256'], 'Source hash mismatch: '+name)

def kernels(root):
    names=['CIE_xyz_1931_2deg.csv','CIE_std_illum_A_1nm.csv',
           'CIE_std_illum_D65.csv','CIE_illum_FLs_1nm.csv']
    tables=[]
    for name in names:
        with (root/'data'/name).open(encoding='utf-8-sig',newline='') as f:
            tables.append({int(row[0]):row for row in csv.reader(f) if row})
    response=[]
    for table,column in [(1,1),(2,1),(3,2),(3,7),(3,11),(3,3),(3,8),(3,12)]:
        light=[]
        for channel in (1,2,3):
            cells=[]
            for wavelength in range(380,780):
                s0,s1=[F(tables[table][w][column]) for w in (wavelength,wavelength+1)]
                f0,f1=[F(tables[0][w][channel]) for w in (wavelength,wavelength+1)]
                require(min(s0,s1,f0,f1)>=0, 'Nonnegative source assumption failed')
                cells.append((s0*f0,(s1-s0)*f0+s0*(f1-f0),(s1-s0)*(f1-f0)))
            light.append(cells)
        normalizer=sum((a+b/2+c/3 for a,b,c in light[1]),F(0))
        require(normalizer>0, 'Nonpositive white normalization')
        response.extend([[(a/normalizer,b/normalizer,c/normalizer) for a,b,c in channel] for channel in light])
    require(all(sum((a+b/2+c/3 for a,b,c in response[j]),F(0))==1 for j in range(1,24,3)), 'White Y identity failed')
    return response

def split_cubic(beta):
    a,b,c,d=beta
    ab,bc,cd=(a+b)/2,(b+c)/2,(c+d)/2
    abc,bcd=(ab+bc)/2,(bc+cd)/2
    mid=(abc+bcd)/2
    return (a,ab,abc,mid),(mid,bcd,cd,d)

def absolute_cubic(beta,width,tolerance):
    """Exact enclosure using positive Bernstein basis and dyadic subdivision."""
    low=high=F(0)
    stack=[(tuple(beta),width,tolerance,0)]
    while stack:
        b,h,tol,depth=stack.pop()
        signed=h*sum(b)/4
        lower=abs(signed)
        upper=h*sum(map(abs,b))/4
        if upper-lower<=tol or depth==48:
            low+=lower;high+=upper
        else:
            left,right=split_cubic(b)
            stack.extend([(right,h/2,tol/2,depth+1),(left,h/2,tol/2,depth+1)])
    require(high-low<=tolerance, 'Cubic integration accuracy not attained')
    return low,high

def shifted(coefficients,offset):
    a,b,c=coefficients
    return a+b*offset+c*offset**2,b+2*c*offset,c

def verify_dual(record,polynomials):
    h=F(record['h']);require(h in (F(1),F(1,2)), 'Unsupported certificate mesh')
    n=int(400/h);v=list(map(F,record['v']));p=list(map(F,record['flux']))
    require(len(v)==15 and len(p)==n+1, 'Dual dimensions')
    require(p[0]==p[-1]==0, 'Dual flux endpoint must be zero')
    q=int(record['coordinate'])-1;require(0<=q<9,'Target coordinate')
    norm=sum(map(abs,v));require(norm==F(record['N']), 'Dual test norm')
    residual=lower=upper=F(0)
    for j in range(n):
        position=j*h;cell=int(position);offset=position-cell
        a,b,c=shifted(polynomials[15+q][cell],offset)
        for i,weight in enumerate(v):
            aa,bb,cc=shifted(polynomials[i][cell],offset)
            a-=weight*aa;b-=weight*bb;c-=weight*cc
        # For u=a+bt+ct^2, construct p with p'=R/h-u and prescribed endpoints.
        integral=a*h+b*h*h/2+c*h**3/3
        R=integral+p[j+1]-p[j]
        residual+=abs(R)
        beta=(p[j],(2*p[j]+p[j+1])/3+b*h*h/6+c*h**3/9,
              (p[j]+2*p[j+1])/3+b*h*h/6+2*c*h**3/9,p[j+1])
        lo,hi=absolute_cubic(beta,h,F(1,10**12)/n)
        lower+=lo;upper+=hi
    require(residual==F(record['R']), 'Dual residual identity failed')
    require(0<upper<=F(record['P']), 'Saved flux upper bound not verified')
    P=F(record['P'])
    require(F(record['a'])==(TAU-residual)/(2*P) and F(record['b'])==-norm/(2*P), 'Dual affine slope formula')
    require(F(record['left'])==0 and F(record['right'])==END, 'Dual parameter domain')
    return dict(id=record['id'],coordinate=q+1,R=residual,P=P,N=norm,a=F(record['a']),b=F(record['b']),left=F(0),right=END)

def nodal_integrals(polynomials,h):
    require(h in (F(1),F(1,2)), 'Unsupported witness mesh')
    n=int(400/h);weights=[]
    for row in polynomials:
        vector=[F(0)]*(n+1)
        for j in range(n):
            position=j*h;cell=int(position)
            a,b,c=shifted(row[cell],position-cell)
            i0=a*h+b*h**2/2+c*h**3/3
            i1=a*h**2/2+b*h**3/3+c*h**4/4
            vector[j]+=i0-i1/h;vector[j+1]+=i1/h
        weights.append(vector)
    return weights

def verify_witness(record,weights):
    h=F(record['h']);d=list(map(F,record['nodes']));n=int(400/h)
    require(h in (F(1),F(1,2)), 'Unsupported witness mesh')
    require(len(d)==n+1, 'Witness node count')
    require(max(map(abs,d))<=1, 'Witness violates passive amplitude')
    g=max(abs(y-x)/(2*h) for x,y in zip(d,d[1:]))
    require(g==F(record['g']), 'Witness slope identity')
    responses=[sum((a*b for a,b in zip(row,d)),F(0)) for row in weights]
    tests=record['tests'];require(len(set(tests))==len(tests) and all(i in range(1,6) for i in tests), 'Witness test IDs')
    selected=[responses[3*(light-1)+c] for light in tests for c in range(3)]
    if record['id'].startswith('R162:'):
        require(max(map(abs,selected))==0, 'Exact four-light matching failed')
        require(g<=F(71,2000), 'Four-light witness outside declared slope class')
        require(max(map(abs,responses[15:]))>TAU, 'Four-light strict failure lost')
    require(selected==list(map(F,record['testResponses'])), 'Witness test integration identity')
    require(responses[15:]==list(map(F,record['targets'])), 'Witness target integration identity')
    return dict(id=record['id'],g=g,delta=max(map(abs,selected)),targets=responses[15:],tests=tests)

def line_value(line,x):return line['a']+line['b']*x

def lines_from_saved(records):
    return [dict(r,**{k:F(r[k]) for k in ('a','b','left','right')}) for r in records]

def verify_failure_lines(proof,witnesses):
    points={}
    for p in proof['points']:
        w=witnesses[p['witnessId']];j=p['coordinate']-1;orientation=p['orientation']
        require(w['tests']==list(range(1,6)), 'Frontier requires all five tests')
        require(orientation in (-1,1), 'Witness sign')
        target=orientation*w['targets'][j]
        require(target>TAU, 'Frontier witness is not a strict failure')
        require(target==F(p['target']) and w['g']==F(p['g']) and w['delta']==F(p['delta']), 'Frontier point identity')
        points[p['id']]=dict(coordinate=j+1,delta=w['delta'],g=w['g'])
    lines=lines_from_saved(proof['failureLines'])
    for line in lines:
        pp=[points[k] for k in line['pointIds']]
        require(len(pp) in (1,2), 'Failure segment must have one or two witnesses')
        require(all(p['coordinate']==line['coordinate'] for p in pp), 'Mixed-coordinate interpolation is invalid')
        if len(pp)==1:
            p=pp[0]
            require(line['a']==p['g'] and line['b']==0 and line['left']==p['delta'] and line['right']==END, 'Constant failure line')
        else:
            p,q=sorted(pp,key=lambda p:p['delta'])
            require(p['delta']<q['delta'], 'Repeated interpolation tolerance')
            require(line_value(line,p['delta'])==p['g'] and line_value(line,q['delta'])==q['g'], 'Witness interpolation identity')
            require(line['left']==max(F(0),p['delta']) and line['right']==min(END,q['delta']), 'Witness interpolation domain')
    return lines,points

def check_partition(pieces):
    require(pieces[0]['left']==0 and pieces[-1]['right']==END, 'Envelope domain coverage')
    for i,p in enumerate(pieces):
        require(p['left']<p['right'], 'Empty envelope piece')
        if i:require(pieces[i-1]['right']==p['left'], 'Gap in envelope partition')

def verify_envelopes(proof,duals,failure_lines):
    safe=lines_from_saved(proof['safePieces']);fail=lines_from_saved(proof['failurePieces'])
    check_partition(safe);check_partition(fail)
    # All crossings of the finite bank suffice to test equality of affine envelopes.
    crossings={F(0),END}
    for x,y in combinations(duals,2):
        if x['b']!=y['b']:
            t=(y['a']-x['a'])/(x['b']-y['b'])
            if 0<t<END:crossings.add(t)
    dual_ids={d['id']:d for d in duals}
    for p in safe:
        source=dual_ids[p['id']]
        require(p['a']==source['a'] and p['b']==source['b'], 'Safe piece is not its cited dual')
        cuts=sorted({p['left'],p['right']}|{x for x in crossings if p['left']<x<p['right']})
        for left,right in zip(cuts,cuts[1:]):
            t=(left+right)/2
            expected=min(max(line_value(d,t) for d in duals if d['coordinate']==j) for j in range(1,10))
            require(line_value(p,t)==expected, 'Safe envelope selection failed')
        require(min(line_value(p,p['left']),line_value(p,p['right']))>0, 'Negative safe slope')
    failure_ids={d['id']:d for d in failure_lines}
    for p in fail:
        source=failure_ids[p['id']]
        require(p['a']==source['a'] and p['b']==source['b'] and source['left']<=p['left']<p['right']<=source['right'], 'Failure piece source/domain')
        for line in failure_lines:
            l=max(p['left'],line['left']);r=min(p['right'],line['right'])
            if l<=r:
                require(all(line_value(p,t)<=line_value(line,t) for t in (l,r)), 'Failure envelope is not minimal')
    knots=sorted({p[k] for p in safe+fail for k in ('left','right')})
    gap=F(0);locations=[]
    for l,r in zip(knots,knots[1:]):
        t=(l+r)/2
        s=next(p for p in safe if p['left']<=t<=p['right'])
        f=next(p for p in fail if p['left']<=t<=p['right'])
        for x in (l,r):
            width=line_value(f,x)-line_value(s,x)
            require(width>0, 'Safe and failure certificates conflict')
            if width>gap:gap=width;locations=[x]
            elif width==gap:locations.append(x)
    require(gap==F(proof['maximumGap']) and F(proof['maximumGapDelta']) in locations, 'Maximum continuum gap identity')
    return dict(safePieces=len(safe),failurePieces=len(fail),maximumGap=float(gap),uniformGapTargetMet=gap<=F(1,1000))

def verify_count(proof,duals,witnesses):
    claim=proof['minimumCountClaim'];g=F(claim['g']);delta=F(claim['safeDelta'])
    require(g==F(71,2000) and delta==F('0.000181293988') and claim['count']==5, 'Unexpected minimum-count claim')
    upper=max(min(d['R']+2*g*d['P']+delta*d['N'] for d in duals if d['coordinate']==j) for j in range(1,10))
    require(upper<=TAU, 'Five-light positive-tolerance sufficiency failed')
    exclusions=[witnesses[k] for k in proof['excludedFourLightWitnessIds']]
    require(len(exclusions)==5, 'Four-light exclusion count')
    covered=0
    for n in range(5):
        for protocol in combinations(range(1,6),n):
            matches=[w for w in exclusions if set(protocol)<=set(w['tests']) and w['g']<=g and w['delta']==0 and max(map(abs,w['targets']))>TAU]
            require(bool(matches), 'Uncovered smaller protocol: '+str(protocol));covered+=1
    require(covered==31, 'Finite protocol enumeration')
    return dict(minimumConstituentLights=5,excludedProtocols=covered,g=float(g),pairwiseToleranceThrough=float(delta),fiveLightUpper=float(upper))

def verify(root,progress=True):
    start=time.monotonic();proof=read(root/'proof.json');check_sources(root,proof)
    polynomial=kernels(root);duals=[]
    for i,c in enumerate(proof['duals'],1):
        duals.append(verify_dual(c,polynomial))
        if progress and i%9==0:print(f'Checked {i}/{len(proof["duals"])} continuous dual certificates',flush=True)
    weights={F(h):nodal_integrals(polynomial,F(h)) for h in {r['h'] for r in proof['witnesses']}}
    witnesses={r['id']:verify_witness(r,weights[F(r['h'])]) for r in proof['witnesses']}
    if progress:print(f'Checked {len(witnesses)} exact continuous witnesses',flush=True)
    failure_lines,points=verify_failure_lines(proof,witnesses)
    frontier=verify_envelopes(proof,duals,failure_lines);count=verify_count(proof,duals,witnesses)
    return dict(status='PASS',checkedAtUtc=datetime.now(timezone.utc).isoformat(),optimizerCalls=0,arithmetic='Exact rational identities and rigorous cubic integral enclosures',dualCount=len(duals),witnessCount=len(witnesses),strictFrontierPoints=len(points),frontier=frontier,testCount=count,proofSHA256=digest(root/'proof.json'),elapsedSeconds=round(time.monotonic()-start,3))

if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--output',type=Path);args=parser.parse_args()
    result=verify(Path(__file__).resolve().parent)
    if args.output:args.output.write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8')
    print(json.dumps(result,indent=2))
