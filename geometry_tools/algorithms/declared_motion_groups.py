"""Close explicit rigid motion dependencies without welding source geometry."""
import numpy as np

def declared_motion_groups(maps,links,*,fixed_groups=()):
    maps=[np.asarray(x,int) for x in maps]
    count=max(int(x.max()) for x in maps)+1
    parent=np.arange(count)
    def root(i):
        while parent[i]!=i:parent[i]=parent[parent[i]];i=parent[i]
        return int(i)
    closures=[]
    for link in links:
        groups=np.unique(np.r_[maps[link['firstActor']][link['firstControls']],maps[link['secondActor']][link['secondControls']]])
        assert len(groups)>1 and groups.min()>=0
        seed=root(int(groups[0]))
        for g in groups[1:]:parent[root(int(g))]=seed
        closures.append(dict(sourceGroups=groups.tolist(),reason=link['reason']))
    roots=np.array([root(i) for i in range(count)])
    fixed=np.asarray(list(fixed_groups),int)
    for closure in closures:
        if len(fixed) and np.any(roots[fixed]==roots[closure['sourceGroups'][0]]):
            raise ValueError('Declared moving interface overlaps a separately fixed source group')
    _,normalized=np.unique(roots,return_inverse=True)
    result=[normalized[x] for x in maps]
    return result,normalized,dict(complete=True,originalGroups=count,mergedGroups=int(normalized.max())+1,
        explicitClosures=closures,positionsModified=False,sourcePointsWelded=False,
        rule='One shared translation per explicitly declared source interface star; original point offsets are retained',
        anatomicalAttachmentInferred=False)
