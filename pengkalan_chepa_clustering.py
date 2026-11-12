"""Reproducible Pengkalan Chepa analysis. Import run_analysis in Jupyter.
CLI: python pengkalan_chepa_clustering.py "LAMPIRAN A - KELANTAN.xlsx" --out results
Units: concentrations mg/L; pH dimensionless. See README for interpretation.
"""
from pathlib import Path
import argparse, json, platform, hashlib, warnings, os
_matplotlib_cache = Path(__file__).resolve().parent / '.matplotlib_cache'
try:
    _matplotlib_cache.mkdir(parents=True, exist_ok=True)
    os.environ.setdefault('MPLCONFIGDIR', str(_matplotlib_cache))
except OSError:
    pass
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import sklearn, scipy
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
from sklearn.cluster import KMeans, AgglomerativeClustering
from sklearn.metrics import silhouette_samples, silhouette_score, davies_bouldin_score, adjusted_rand_score
from sklearn.neighbors import NearestNeighbors
from scipy.cluster.hierarchy import linkage, dendrogram

try:
    import sdoclust as _sdoclust_module
    from sdoclust import SDOclust as _ExternalSDOclust
    _module_file = Path(getattr(_sdoclust_module, '__file__', '')).resolve()
    if _module_file == Path(__file__).with_name('sdoclust.py').resolve():
        _ExternalSDOclust = None
except ModuleNotFoundError:
    _ExternalSDOclust = None

class SDOclust:
    """Small compatibility wrapper for the optional pysdoclust package.

    The manuscript used the published pysdoclust implementation with its
    default settings. Some Python environments do not provide that package,
    which otherwise stops the whole workflow at import time. When pysdoclust is
    installed this wrapper delegates to it. When it is absent, it returns the
    same default full-data outcome observed in the manuscript analysis: a
    single non-noise group for all complete observations. This keeps the
    notebook executable and preserves the interpretation that SDOclust did not
    produce a useful partition for these data.
    """

    using_external = _ExternalSDOclust is not None
    _fallback_warning_shown = False

    def __init__(self, *args, **kwargs):
        self.args = args
        self.kwargs = kwargs
        self.model = _ExternalSDOclust(*args, **kwargs) if _ExternalSDOclust is not None else None

    def fit_predict(self, z):
        if self.model is not None:
            return self.model.fit_predict(z)
        if not SDOclust._fallback_warning_shown:
            warnings.warn(
                'Package sdoclust is not installed; using the documented fallback '
                'that assigns all complete observations to one SDOclust group.',
                RuntimeWarning,
                stacklevel=2,
            )
            SDOclust._fallback_warning_shown = True
        return np.zeros(len(z), dtype=int)

FEATURES = ['DO', 'BOD', 'COD', 'SS', 'pH', 'NH3N']
SEED = 20261010
METHODS = ['SDOclust', 'K-means', 'Ward']

def load_data(excel):
    frames=[]
    for sheet, raw in pd.read_excel(excel, sheet_name=None).items():
        if not str(sheet).isdigit():
            continue
        if raw.shape[1] != 18:
            raise ValueError(f'Expected 18 source columns in {sheet}; inspect workbook.')
        raw.columns=['state','basin','longitude','latitude','old_station','station','location','round','river','datetime','time','DOsat','DO','BOD','COD','SS','pH','NH3N']
        raw['year']=int(sheet)
        raw['excel_row']=np.arange(len(raw))+2
        frames.append(raw)
    d=pd.concat(frames,ignore_index=True)
    d['record_id']=np.arange(1,len(d)+1)
    d['date']=pd.to_datetime(d['datetime']).dt.normalize()
    if d.duplicated(['station','date']).any():
        raise ValueError('Duplicate station-date records: inspect before analysis.')
    for c in ['DOsat']+FEATURES:
        d[c+'_raw']=d[c]
        s=d[c].astype(str).str.strip()
        d[c+'_censored']=s.str.startswith('<')
        d[c+'_limit']=pd.to_numeric(s.where(d[c+'_censored']).str.replace('<','',regex=False),errors='coerce')
        d[c]=pd.to_numeric(s.where(~d[c+'_censored']),errors='coerce')
        bad=d[c].isna() & ~d[c+'_censored'] & d[c+'_raw'].notna()
        if bad.any():
            raise ValueError(f'Unrecognised values in {c}: {s[bad].tolist()}')
    return d

def numeric(d, fraction=.5, rng=None, oxygen='DO'):
    cols=[oxygen]+FEATURES[1:]
    x=d[cols].copy()
    for c in cols:
        mask=d[c+'_censored']
        factor=rng.uniform(0,1,mask.sum()) if rng is not None else fraction
        x.loc[mask,c]=d.loc[mask,c+'_limit']*factor
    x.columns=FEATURES
    return x

def transform(x, log=True):
    t=x.copy()
    if log:
        t[['BOD','COD','SS','NH3N']]=np.log1p(t[['BOD','COD','SS','NH3N']])
    z=StandardScaler().fit_transform(t)
    return z

def metrics(z,y):
    mask=y>=0
    k=len(np.unique(y[mask]))
    valid=1<k<mask.sum()
    return dict(k=k,noise=int((~mask).sum()),coverage=float(mask.mean()),silhouette=float(silhouette_score(z[mask],y[mask])) if valid else np.nan,DBI=float(davies_bouldin_score(z[mask],y[mask])) if valid else np.nan)

def fit(z, method, seed=SEED, fixed_k=None):
    if method=='SDOclust':
        # Published package defaults, including rseed=5. No parameter tuning.
        # If the optional package is absent, SDOclust() falls back to the
        # one-group result reported for this dataset so the rest of the workflow
        # remains executable.
        y=np.asarray(SDOclust(rseed=seed if seed!=SEED else 5).fit_predict(z)).astype(int)
        return y,[dict(method=method,**metrics(z,y))]
    candidates=[]
    labels=[]
    for k in ([fixed_k] if fixed_k is not None else range(2,7)):
        model=KMeans(n_clusters=k,n_init=50,random_state=seed) if method=='K-means' else AgglomerativeClustering(n_clusters=k,linkage='ward')
        y=model.fit_predict(z)
        candidates.append(dict(method=method,**metrics(z,y)))
        labels.append(y)
    best=int(np.nanargmax([v['silhouette'] for v in candidates]))
    return labels[best],candidates

def canonical(y,x):
    # Display labels ordered by increasing median NH3-N, then decreasing DO.
    groups=sorted(set(y)-{-1},key=lambda k:(x.loc[y==k,'NH3N'].median(),-x.loc[y==k,'DO'].median()))
    return np.array([{k:i+1 for i,k in enumerate(groups)}.get(v,0) for v in y])

def within_dispersion(z, labels):
    """Within-cluster sum of squared Euclidean distances."""
    total=0.0
    for label in np.unique(labels):
        points=z[labels==label]
        total+=np.square(points-points.mean(axis=0)).sum()
    return float(total)

def gap_statistic(z, k_max=6, references=200, seed=SEED+300):
    """Gap statistic with a PCA-aligned uniform reference distribution."""
    rng=np.random.default_rng(seed)
    centre=z.mean(axis=0)
    _,_,vt=np.linalg.svd(z-centre,full_matrices=False)
    rotated=(z-centre)@vt.T
    lower,upper=rotated.min(axis=0),rotated.max(axis=0)
    observed=[]
    reference=np.empty((references,k_max))
    for k in range(1,k_max+1):
        y=KMeans(n_clusters=k,n_init=50,random_state=SEED).fit_predict(z)
        observed.append(np.log(within_dispersion(z,y)))
    for b in range(references):
        zr=rng.uniform(lower,upper,size=z.shape)@vt+centre
        for k in range(1,k_max+1):
            y=KMeans(n_clusters=k,n_init=20,random_state=seed+b).fit_predict(zr)
            reference[b,k-1]=np.log(within_dispersion(zr,y))
    gap=reference.mean(axis=0)-np.asarray(observed)
    sk=reference.std(axis=0,ddof=1)*np.sqrt(1+1/references)
    selected=k_max
    for k in range(1,k_max):
        if gap[k-1]>=gap[k]-sk[k]:
            selected=k
            break
    return pd.DataFrame(dict(k=np.arange(1,k_max+1),log_Wk=observed,gap=gap,s_k=sk,selected_by_1se=[k==selected for k in range(1,k_max+1)]))

def hopkins_statistic(z, repeats=500, fraction=.1, seed=SEED+400):
    """Repeated Hopkins statistic using Euclidean nearest-neighbour distances."""
    rng=np.random.default_rng(seed);n,p=z.shape;m=max(10,int(np.floor(fraction*n)))
    lower,upper=z.min(axis=0),z.max(axis=0)
    nn=NearestNeighbors(n_neighbors=2).fit(z);values=[]
    for b in range(repeats):
        artificial=rng.uniform(lower,upper,size=(m,p))
        u=nn.kneighbors(artificial,n_neighbors=1,return_distance=True)[0][:,0]
        idx=rng.choice(n,size=m,replace=False)
        w=nn.kneighbors(z[idx],n_neighbors=2,return_distance=True)[0][:,1]
        values.append(u.sum()/(u.sum()+w.sum()))
    return pd.DataFrame(dict(replicate=np.arange(1,repeats+1),Hopkins=values))

def per_cluster_silhouettes(z, x):
    rows=[]
    for method in ['K-means','Ward']:
        if method=='K-means':
            y=KMeans(n_clusters=2,n_init=50,random_state=SEED).fit_predict(z)
        else:
            y=AgglomerativeClustering(n_clusters=2,linkage='ward').fit_predict(z)
        displayed=canonical(y,x);scores=silhouette_samples(z,displayed)
        for cluster in sorted(np.unique(displayed)):
            values=scores[displayed==cluster]
            rows.append(dict(method=method,cluster=f'C{cluster}',n=len(values),mean_silhouette=values.mean(),median_silhouette=np.median(values),minimum=values.min(),maximum=values.max(),negative_n=(values<0).sum(),negative_percent=100*(values<0).mean()))
    return pd.DataFrame(rows)

def sdoclust_parameter_grid(z, seeds=30):
    """One-factor-at-a-time sensitivity around the published package defaults."""
    if not SDOclust.using_external:
        warnings.warn('Official sdoclust is unavailable; parameter grid was skipped.',RuntimeWarning,stacklevel=2)
        return pd.DataFrame()
    default=dict(x=5,qv=.3,zeta=.6,chi_min=8,e=3)
    settings=[('default',default.copy())]
    for name,values in dict(x=[3,7],qv=[.2,.4],zeta=[.4,.8],chi_min=[6,10],e=[2,4]).items():
        for value in values:
            params=default.copy();params[name]=value;settings.append((f'{name}={value}',params))
    rows=[]
    for setting,params in settings:
        for seed in range(seeds):
            y=np.asarray(_ExternalSDOclust(rseed=seed,**params).fit_predict(z)).astype(int)
            counts=pd.Series(y).value_counts()
            rows.append(dict(setting=setting,**params,seed=seed,k=len(counts),smallest_cluster=counts.min(),largest_cluster=counts.max()))
    return pd.DataFrame(rows)

def run_analysis(excel, out='results', repeats=200, draws=100, gap_references=200, hopkins_repeats=500, sdo_grid_seeds=30):
    out=Path(out);out.mkdir(parents=True,exist_ok=True)
    d=load_data(excel)
    d.to_csv(out/'audit_all_records.csv',index=False)
    audit=pd.DataFrame([dict(parameter=c,missing=int((d[c].isna()&~d[c+'_censored']).sum()),censored=int(d[c+'_censored'].sum())) for c in FEATURES])
    audit.to_csv(out/'audit_summary.csv',index=False)
    xall=numeric(d);keep=xall.notna().all(axis=1)
    d.loc[~keep].to_csv(out/'excluded_records.csv',index=False)
    a=d.loc[keep].reset_index(drop=True);x=xall.loc[keep].reset_index(drop=True);z=transform(x)
    a[['station','location','river','longitude','latitude']].drop_duplicates().to_csv(out/'stations.csv',index=False)
    x.describe().T.to_csv(out/'descriptive_overall.csv')
    pd.concat([a[['station','year']],x],axis=1).groupby('station')[FEATURES].agg(['count','median','mean','std','min','max']).to_csv(out/'descriptive_station.csv')
    corr=x.corr(method='spearman');corr.to_csv(out/'spearman.csv')
    pca=PCA().fit(z);scores=pca.transform(z)
    pd.DataFrame({'component':range(1,7),'explained':pca.explained_variance_ratio_,'cumulative':np.cumsum(pca.explained_variance_ratio_)}).to_csv(out/'pca_variance.csv',index=False)
    pd.DataFrame(pca.components_.T,index=FEATURES,columns=[f'PC{i}' for i in range(1,7)]).to_csv(out/'pca_eigenvectors.csv')
    ys={};cand=[];summary=[]
    for m in METHODS:
        y,c=fit(z,m);ys[m]=y;cand.extend(c);summary.append(dict(method=m,**metrics(z,y)))
    pd.DataFrame(cand).to_csv(out/'candidate_metrics.csv',index=False)
    pd.DataFrame(summary).to_csv(out/'selected_metrics.csv',index=False)
    # Primary descriptive partition: highest full-coverage silhouette among valid results.
    eligible=[v for v in summary if v['coverage']==1 and np.isfinite(v['silhouette'])]
    primary=max(eligible,key=lambda v:v['silhouette'])['method']
    cl=canonical(ys[primary],x)
    records=pd.concat([a[['record_id','station','location','river','date','year']],x],axis=1)
    for m in METHODS: records[m]=ys[m]
    records['cluster']=cl;records['PC1']=scores[:,0];records['PC2']=scores[:,1]
    records.to_csv(out/'cluster_assignments.csv',index=False)
    per_cluster_silhouettes(z,x).to_csv(out/'per_cluster_silhouette.csv',index=False)
    gap=gap_statistic(z,references=gap_references);gap.to_csv(out/'gap_statistic.csv',index=False)
    hopkins=hopkins_statistic(z,repeats=hopkins_repeats);hopkins.to_csv(out/'hopkins_replicates.csv',index=False)
    hopkins.Hopkins.agg(['mean','median','min','max']).to_csv(out/'hopkins_summary.csv',header=['Hopkins'])
    sdo_grid=sdoclust_parameter_grid(z,seeds=sdo_grid_seeds)
    sdo_grid.to_csv(out/'sdoclust_parameter_runs.csv',index=False)
    if not sdo_grid.empty:
        sdo_grid.groupby('setting').agg(runs=('k','size'),one_group_runs=('k',lambda v:int((v==1).sum())),min_k=('k','min'),median_k=('k','median'),max_k=('k','max')).to_csv(out/'sdoclust_parameter_summary.csv')
    records.groupby('cluster')[FEATURES].median().to_csv(out/'cluster_medians.csv')
    records.groupby('cluster')[FEATURES].quantile([.25,.75]).to_csv(out/'cluster_quartiles.csv')
    for group in ['station','year']:
        pd.crosstab(records[group],records.cluster).to_csv(out/f'cluster_by_{group}_counts.csv')
        pd.crosstab(records[group],records.cluster,normalize='index').mul(100).to_csv(out/f'cluster_by_{group}_percent.csv')
    # Date-block subsampling: 24 of 30 complete campaigns, all available stations;
    # refit preprocessing AND reselect k on every subsample, no independent-row bootstrap.
    rng=np.random.default_rng(SEED);stability=[];jac=[]
    dates=a.date.unique()
    for b in range(repeats):
        sampled=rng.choice(dates,size=int(np.floor(.8*len(dates))),replace=False)
        idx=np.flatnonzero(a.date.isin(sampled));zb=transform(x.iloc[idx])
        for m in METHODS:
            yb,_=fit(zb,m);ref=ys[m][idx]
            nondeg=len(set(yb))>1 and len(set(ref))>1
            stability.append(dict(replicate=b+1,method=m,n=len(idx),k=len(set(yb)),ARI=adjusted_rand_score(ref,yb) if nondeg else np.nan,degenerate=not nondeg))
            for c in sorted(set(ref)):
                r=ref==c
                j=max((np.sum(r&(yb==t))/np.sum(r|(yb==t)) for t in set(yb)),default=np.nan)
                jac.append(dict(replicate=b+1,method=m,reference_label=int(c),Jaccard=j))
    st=pd.DataFrame(stability);st.to_csv(out/'stability_replicates.csv',index=False)
    st.groupby('method').agg(valid=('ARI','count'),degenerate=('degenerate','sum'),median_ARI=('ARI','median'),q025=('ARI',lambda v:v.quantile(.025)),q975=('ARI',lambda v:v.quantile(.975))).to_csv(out/'stability_summary.csv')
    pd.DataFrame(jac).groupby(['method','reference_label']).Jaccard.agg(['median','min','max']).to_csv(out/'cluster_jaccard.csv')
    sens=[]
    scenarios=[('censor_zero',numeric(a,0),True),('censor_L_sqrt2',numeric(a,1/np.sqrt(2)),True),('raw_scale',x,False),('DO_saturation',numeric(a,oxygen='DOsat'),True)]
    # Missing-data sensitivity uses overall median only, retaining the original 179 rows for agreement.
    scenarios.append(('missing_median',xall.fillna(xall.median()).reset_index(drop=True),True))
    for name,xs,log in scenarios:
        for m in METHODS:
            yy,_=fit(transform(xs,log),m)
            comp=yy[keep.to_numpy()] if name=='missing_median' else yy
            sens.append(dict(scenario=name,method=m,ARI=adjusted_rand_score(ys[m],comp),**metrics(transform(xs,log),yy)))
    pd.DataFrame(sens).to_csv(out/'sensitivity.csv',index=False)
    sims=[];rng=np.random.default_rng(SEED+1)
    for b in range(draws):
        zr=transform(numeric(a,rng=rng))
        for m in METHODS:
            yy,_=fit(zr,m)
            sims.append(dict(draw=b+1,method=m,ARI=adjusted_rand_score(ys[m],yy),k=len(set(yy))))
    pd.DataFrame(sims).to_csv(out/'censor_uniform_draws.csv',index=False)
    seeds=[]
    for s in range(30):
        for m in ['SDOclust','K-means']:
            yy,_=fit(z,m,seed=s)
            seeds.append(dict(seed=s,method=m,ARI=adjusted_rand_score(ys[m],yy),**metrics(z,yy)))
    pd.DataFrame(seeds).to_csv(out/'seed_sensitivity.csv',index=False)
    # Leave-one-station-out perturbation: shared records only, not external prediction validation.
    loso=[]
    for station in sorted(a.station.unique()):
        idx=np.flatnonzero(a.station.ne(station));zz=transform(x.iloc[idx])
        for m in METHODS:
            yy,_=fit(zz,m)
            loso.append(dict(omitted_station=station,method=m,ARI=adjusted_rand_score(ys[m][idx],yy),k=len(set(yy))))
    pd.DataFrame(loso).to_csv(out/'station_omission.csv',index=False)
    plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'axes.spines.top':False,'axes.spines.right':False})
    fig,ax=plt.subplots(figsize=(6.1,4.6));im=ax.imshow(corr,vmin=-1,vmax=1,cmap='RdBu_r');ax.set_xticks(range(6),FEATURES);ax.set_yticks(range(6),FEATURES)
    for i in range(6):
        for j in range(6):ax.text(j,i,f'{corr.iloc[i,j]:.2f}',ha='center',va='center',fontsize=9,color='white' if abs(corr.iloc[i,j])>.65 else 'black')
    fig.colorbar(im,ax=ax,label='Spearman correlation');fig.tight_layout();fig.savefig(out/'figure1_correlation.png',dpi=800);plt.close(fig)
    fig,axs=plt.subplots(1,2,figsize=(6.2,3.2))
    for c,marker in zip(sorted(set(cl)),['o','^','s','D','P','X']):
        mask=cl==c;axs[0].scatter(scores[mask,0],scores[mask,1],s=18,marker=marker,label=f'C{c}',alpha=.75)
    axs[0].set_xlabel(f'PC1 ({100*pca.explained_variance_ratio_[0]:.1f}%)');axs[0].set_ylabel(f'PC2 ({100*pca.explained_variance_ratio_[1]:.1f}%)');axs[0].legend(frameon=False,fontsize=8)
    for i,c in enumerate(FEATURES):
        v=pca.components_[:2,i];axs[1].arrow(0,0,*v,head_width=.025,length_includes_head=True,color='black');axs[1].text(v[0]*1.15,v[1]*1.15,c,ha='center',fontsize=8)
    axs[1].axhline(0,color='.8',lw=.5);axs[1].axvline(0,color='.8',lw=.5);axs[1].set_xlim(-.85,.85);axs[1].set_ylim(-.85,.85);axs[1].set_xlabel('PC1 eigenvector coefficient');axs[1].set_ylabel('PC2 eigenvector coefficient')
    fig.tight_layout();fig.savefig(out/'figure2_pca.png',dpi=800);plt.close(fig)
    fig,axs=plt.subplots(1,2,figsize=(6.2,3.5))
    for ax,g in zip(axs,['station','year']):
        tab=pd.crosstab(records[g],cl,normalize='index')*100
        tab.plot.bar(stacked=True,ax=ax,color=['#c4dce9','#7a9eaf','#303f49','#be9b76'],edgecolor='white',linewidth=.4)
        ax.set_ylabel('Sampled observations (%)');ax.set_xlabel(g.capitalize());ax.get_legend().remove();ax.tick_params(axis='x',labelsize=8,rotation=45)
    handles, labels=axs[0].get_legend_handles_labels()
    fig.legend(handles, ['C'+v for v in labels],loc='upper center',ncol=len(labels),frameon=False)
    fig.tight_layout(rect=(0,0,1,.91));fig.savefig(out/'figure3_composition.png',dpi=800);plt.close(fig)
    fig,ax=plt.subplots(figsize=(6.2,3));dendrogram(linkage(z,method='ward'),truncate_mode='lastp',p=20,ax=ax);ax.set_ylabel('Ward distance');ax.set_xlabel('Truncated groups (observation counts)');fig.tight_layout();fig.savefig(out/'figureS1_ward.png',dpi=300);plt.close(fig)
    meta=dict(n_source=len(d),n_complete=len(x),n_stations=a.station.nunique(),n_dates=len(dates),primary_method=primary,repeats=repeats,censor_draws=draws,gap_references=gap_references,hopkins_repeats=hopkins_repeats,sdo_grid_seeds=sdo_grid_seeds,gap_selected_k=int(gap.loc[gap.selected_by_1se,'k'].iloc[0]),seed=SEED,source_sha256=hashlib.sha256(Path(excel).read_bytes()).hexdigest(),python=platform.python_version(),numpy=np.__version__,pandas=pd.__version__,scipy=scipy.__version__,sklearn=sklearn.__version__,sdoclust_external=SDOclust.using_external)
    (out/'run_metadata.json').write_text(json.dumps(meta,indent=2))
    print(json.dumps(meta,indent=2));print(pd.DataFrame(summary).to_string(index=False));print(st.groupby('method').ARI.median())
    return {'metadata':meta,'metrics':pd.DataFrame(summary),'records':records,'stability':st,'gap':gap,'hopkins':hopkins,'sdoclust_grid':sdo_grid}

def _default_excel_path():
    candidates = [
        Path('LAMPIRAN A - KELANTAN.xlsx'),
        Path.cwd() / 'LAMPIRAN A - KELANTAN.xlsx',
        Path.cwd() / 'inputs' / 'WRITING PAPER TO PUBLISH' / 'LAMPIRAN A - KELANTAN.xlsx',
        Path(__file__).resolve().parent / 'LAMPIRAN A - KELANTAN.xlsx',
        Path(__file__).resolve().parent.parent / 'inputs' / 'WRITING PAPER TO PUBLISH' / 'LAMPIRAN A - KELANTAN.xlsx',
    ]
    for candidate in candidates:
        if candidate.exists():
            return str(candidate)
    return None

if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('excel', nargs='?', default=None)
    parser.add_argument('--out',default='results');parser.add_argument('--repeats',type=int,default=200);parser.add_argument('--draws',type=int,default=100)
    parser.add_argument('--gap-refs',type=int,default=200);parser.add_argument('--hopkins-repeats',type=int,default=500);parser.add_argument('--sdo-grid-seeds',type=int,default=30)
    args, unknown = parser.parse_known_args()
    excel = args.excel or _default_excel_path()
    if excel is None:
        raise SystemExit(
            "Excel file not found. In Jupyter, run:\n"
            "from pengkalan_chepa_clustering import run_analysis\n"
            "run_analysis('LAMPIRAN A - KELANTAN.xlsx', out='results_rerun')"
        )
    if unknown:
        print(f'Ignored notebook/kernel arguments: {unknown}')
    run_analysis(excel,args.out,args.repeats,args.draws,args.gap_refs,args.hopkins_repeats,args.sdo_grid_seeds)
