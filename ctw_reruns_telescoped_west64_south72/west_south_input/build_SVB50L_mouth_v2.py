from pathlib import Path
import numpy as np
from scipy.ndimage import label
import matplotlib;matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
R=Path('/mnt/data');ny,nx=684,576
read=lambda fn:np.fromfile(R/fn,dtype='>f8').reshape(ny,nx)
B=read('bathy_realistic_minusSVB_W64S72.bin'); A=read('bathy_realistic_plusSVB_W64S72.bin')
dx=np.fromfile(R/'delX_W64S72.bin',dtype='>f8');dy=np.fromfile(R/'delY_W64S72.bin',dtype='>f8')
import re
T=(R/'data').read_text();get=lambda key:float(re.search(r'\b'+key+r'\s*=\s*([-+\d.]+)',T).group(1))
x=get('xgOrigin')+np.cumsum(dx)-dx/2;y=get('ygOrigin')+np.cumsum(dy)-dy/2
coast=lambda H: np.array([np.where(row<0)[0].max() for row in H]);cb=coast(B);ca=coast(A)
active=np.where(ca>cb)[0]; orig_n=int(active.max()-active.min()+1)
# Approx half the original latitudinal mouth: round 141/2 -> 71 rows.
# Center near original indentation centroid, not geographical midpoint of outer difference.
ind=np.maximum(0,ca-cb);center=int(round(np.average(np.arange(ny),weights=ind**2)))
width=int(round(orig_n/2)); start=center-width//2;end=start+width-1
# Rounded broad bay: half-sine raised to exponent < 1 produces a broad, rounded head,
# avoids a triangular cross-shore indentation profile.
maxindent=54  # rounded bay head, safely inside eastern boundary
max_coast_col=nx-1-20   # 14 dry columns before eastern model edge
cc=cb.copy()
for j in range(start,end+1):
    u=(j-start)/(width-1)
    amp=np.sin(np.pi*u)**0.85
    proposed=cb[j]+int(round(maxindent*amp))
    cc[j]=min(proposed,max_coast_col)
cc=np.maximum(cc,cb)
# construct profiles from ORIGINAL SVB, keeping the same NoSVB fixed offshore bathymetry.
N=B.copy()
for j in range(start,end+1):
    if cc[j]<=cb[j]:continue
    u=(j-start)/(width-1)
    src=int(round(active.min()+u*(active.max()-active.min())))
    src=np.clip(src,0,ny-1)
    a=int(np.where(B[j]<=-150)[0].max())
    srca=int(np.where(A[src]<=-150)[0].max()); srcc=int(ca[src])
    if cc[j]<=a+1 or srcc<=srca+1:continue
    ids=np.arange(a+1,cc[j]+1)
    fract=(x[ids]-x[a])/(x[cc[j]]-x[a])
    srcx=x[srca]+fract*(x[srcc]-x[srca])
    prof=np.interp(srcx,x[srca:srcc+1],A[src,srca:srcc+1])
    # For near-boundary profile cells monotonically approach shallow shoreline:
    # blend with an analytic coastal depth only if source has isolated poorly-ordered values.
    prof=np.clip(prof,-149,-2)
    # avoid abrupt transition immediately adjacent to protected isobath
    taper=np.minimum(1,(ids-a)/6)
    taper=taper*taper*(3-2*taper)
    prof=np.clip(-149+(prof+149)*taper,-149,-2)
    allowed=B[j,ids]>-150
    N[j,ids[allowed]]=prof[allowed]
# Explicitly preserve NoSVB protected depths, 150 m and offshore region, and ALL boundary cells
protected=B<=-150
assert np.array_equal(N[protected],B[protected])
assert np.array_equal(N<=-150, protected)
assert np.all(N[:,max_coast_col+1:]==B[:,max_coast_col+1:])
assert np.array_equal((N<0)[:,0],(B<0)[:,0])
assert np.array_equal((N<0)[:,-1],(B<0)[:,-1])
assert np.array_equal((N<0)[0],(B<0)[0])
assert np.array_equal((N<0)[-1],(B<0)[-1])
# Check geography of all-new coastline where NoSVB is land
cn=coast(N); added=np.where(cn>cb)[0]
assert len(added)>0
assert cn.max()<=max(cb.max(),max_coast_col)
labels,ncc=label(N<0)
assert ncc==1, f'{ncc} ocean components'
N.astype('>f8').tofile(R/'bathy_realistic_SVB50L_mouth_W64S72.bin')
# Key QA and plots
pairs=[]
for H,name in ((B,'NoSVB'),(A,'Original SVB'),(N,'Short-mouth SVB')):
    masks=H<0
    pairs.append((name,masks.sum(),coast(H).max()))
print('original rows',active.min(),active.max(), 'span',orig_n)
print('target rows',start,end,'actual coastline change',added.min(),added.max(),'span',len(added))
print('max shore column',cn[added].max(), 'buffer',nx-1-cn[added].max())
print('deep modification count',np.count_nonzero((N!=B)&protected))
print('max abs depth change',np.max(abs(N-B)))
print('components',ncc)
print('properties',pairs)
# proper degrees, crop around bay
jlo=max(0,active.min()-22);jhi=min(ny,active.max()+25);ilo=445;ihi=nx
fig,ax=plt.subplots(figsize=(9,7),constrained_layout=True)
for H,col,ls,labelname in ((B,'#3d4652','-','No bay'),(A,'#2674b8','--','Original SVB'),(N,'#e07628','-','Short-mouth SVB')):
    ax.contour(x[ilo:ihi],y[jlo:jhi],(H[jlo:jhi,ilo:ihi]<0).astype(float),levels=[.5],colors=[col],linestyles=[ls],linewidths=2.1)
ax.axvline(x[-1],color='crimson',lw=1.5,ls=':',label='Eastern grid boundary')
ax.set(xlabel='Longitude (°E)',ylabel='Latitude (°N)',title='Three coastlines — reduced alongshore mouth width')
ax.set_aspect(1/np.cos(np.deg2rad(np.mean(y[jlo:jhi]))))
ax.grid(alpha=.18)
ax.legend(handles=[Line2D([0],[0],color=c,ls=ls,lw=2.2,label=nm) for c,ls,nm in [('#3d4652','-','NoSVB'),('#2674b8','--','Original SVB'),('#e07628','-','SVB short mouth'),('crimson',':','Eastern boundary')]],loc='upper left')
fig.savefig(R/'SVB50L_mouth_three_coastlines.png',dpi=200)
fig,axes=plt.subplots(1,3,figsize=(16,6),sharex=True,sharey=True,constrained_layout=True)
for ax,H,title in zip(axes,[B,A,N],['NoSVB','Original SVB','Short-mouth SVB']):
    pcm=ax.pcolormesh(x[ilo:ihi],y[jlo:jhi],np.ma.masked_where(H[jlo:jhi,ilo:ihi]>=0,-H[jlo:jhi,ilo:ihi]),vmin=0,vmax=250,cmap='Blues',shading='auto')
    ax.contour(x[ilo:ihi],y[jlo:jhi],H[jlo:jhi,ilo:ihi],levels=[-150],colors=['orange'],linewidths=1.8)
    ax.contour(x[ilo:ihi],y[jlo:jhi],(H[jlo:jhi,ilo:ihi]<0).astype(float),levels=[.5],colors=['black'],linewidths=1)
    ax.axvline(x[-1],color='crimson',ls=':',lw=1)
    ax.set(title=title,xlabel='Longitude (°E)')
axes[0].set_ylabel('Latitude (°N)')
fig.colorbar(pcm,ax=axes,label='Water depth (m)',shrink=.75)
fig.savefig(R/'SVB50L_mouth_bathymetry.png',dpi=190)
# gradients compare (inclusive shoreline neighbors excluded)
for H,name in ((A,'SVB'),(N,'new')):
    vals=[]
    for h0,h1 in ((H[start:end],H[start+1:end+1]),(H[start:end+1,:-1],H[start:end+1,1:])):
      wet=(h0<0)&(h1<0)&(h0>-150)&(h1>-150)
      r=np.abs(h1[wet]-h0[wet])/(-h1[wet]-h0[wet]);vals.extend(r)
    print(name,'shallow pair r p99/max',np.percentile(vals,99),max(vals))
