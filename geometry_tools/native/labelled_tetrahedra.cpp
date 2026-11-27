// Native implementation of algorithms/labelled_tetrahedra.py. Calculation
// tetrahedra never become runtime meshes. Shared simplex/label keys retain
// reciprocal coordinates and triangle orientation across all cells.
#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <limits>
#include <map>
#include <set>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <vector>
using ID=uint64_t; using V=std::array<double,3>; using F=std::array<ID,3>;
V add(V a,V b){return {a[0]+b[0],a[1]+b[1],a[2]+b[2]};}
V sub(V a,V b){return {a[0]-b[0],a[1]-b[1],a[2]-b[2]};}
V mul(V a,double b){return {a[0]*b,a[1]*b,a[2]*b};}
double dot(V a,V b){return a[0]*b[0]+a[1]*b[1]+a[2]*b[2];}
V cross(V a,V b){return {a[1]*b[2]-a[2]*b[1],a[2]*b[0]-a[0]*b[2],a[0]*b[1]-a[1]*b[0]};}
V unit(V a){double n=std::sqrt(dot(a,a));if(!(n>0))throw std::runtime_error("zero calculation normal");return mul(a,1/n);}
template<class T>T read(std::ifstream& in){T x;in.read(reinterpret_cast<char*>(&x),sizeof x);if(!in)throw std::runtime_error("truncated calculation input");return x;}
template<class T>void write(std::ofstream& out,const T& x){out.write(reinterpret_cast<const char*>(&x),sizeof x);}
struct Point{V p;std::vector<std::pair<ID,double>> weights;ID ties=0;};
struct CutKey{std::array<ID,4> nodes{};ID labels=0;uint8_t size=0;bool operator==(const CutKey& b)const{return size==b.size&&labels==b.labels&&nodes==b.nodes;}};
struct CutHash{size_t operator()(const CutKey& k)const{size_t h=k.labels^k.size;for(auto n:k.nodes)h^=n+0x9e3779b97f4a7c15ULL+(h<<6)+(h>>2);return h;}};
struct FaceHash{size_t operator()(const F& k)const{size_t h=0;for(auto n:k)h^=n+0x9e3779b97f4a7c15ULL+(h<<6)+(h>>2);return h;}};
// Reorthogonalized thin QR, at most four columns. Keeping score equations
// separate from coordinates avoids an arbitrary world-space weld tolerance.
std::vector<double> least_squares(const std::vector<std::vector<double>>& A,const std::vector<double>& b){
 size_t m=A.size(),n=A[0].size();if(m<n)throw std::runtime_error("underdetermined tied score field");
 std::vector<std::vector<double>> Q(n,std::vector<double>(m)),R(n,std::vector<double>(n));
 for(size_t j=0;j<n;j++){
  for(size_t k=0;k<m;k++)Q[j][k]=A[k][j];
  for(int pass=0;pass<2;pass++)for(size_t i=0;i<j;i++){
   double r=0;for(size_t k=0;k<m;k++)r+=Q[i][k]*Q[j][k];R[i][j]+=r;for(size_t k=0;k<m;k++)Q[j][k]-=r*Q[i][k];
  }
  double r=0;for(double x:Q[j])r+=x*x;r=std::sqrt(r);if(!(r>1e-18))throw std::runtime_error("rank-deficient tied score field");R[j][j]=r;for(double& x:Q[j])x/=r;
 }
 std::vector<double> x(n);for(size_t j=0;j<n;j++)for(size_t k=0;k<m;k++)x[j]+=Q[j][k]*b[k];
 for(int j=int(n)-1;j>=0;j--){for(size_t k=j+1;k<n;k++)x[j]-=R[j][k]*x[k];x[j]/=R[j][j];}return x;
}
struct Engine{
 size_t labels;std::vector<Point> points;std::vector<std::vector<double>> score;
 size_t planes=0;std::vector<int64_t> zero_side;std::vector<ID> cell_domain,surface_refinements;std::vector<std::vector<std::pair<ID,int64_t>>> constraints;
 std::unordered_map<CutKey,ID,CutHash> cache;std::map<std::vector<ID>,ID> centroids;
 std::vector<std::unordered_map<F,F,FaceHash>> faces;double max_residual=0,max_weight_error=0,max_pruned_weight=0;ID omitted_zero_area=0;
 double value(ID p,size_t a,size_t b)const{
  const Point& q=points[p];if((q.ties&(ID(1)<<a))&&(q.ties&(ID(1)<<b)))return 0;
  double d=0;for(auto [i,w]:q.weights)d+=w*(score[i][a]-score[i][b]);return d;
 }
 double plane_value(ID p,size_t plane)const{
  const Point& q=points[p];if(q.ties&(ID(1)<<plane))return 0;
  // Two fields may coincide only on this minimal source face/edge. A cut
  // already on one then lies exactly on the other, even if evaluating its
  // solved barycentric weights would leave a tiny arithmetic residual.
  for(size_t tied=0;tied<planes;tied++)if(q.ties&(ID(1)<<tied)){
   bool identical=true;for(auto [i,w]:q.weights)if(score[i][plane]!=score[i][tied]){identical=false;break;}
   if(identical)return 0;
  }
  double d=0,magnitude=0;for(auto [i,w]:q.weights){double term=w*score[i][plane];d+=term;magnitude+=std::abs(term);}
  // An affine identity evaluated through a solved barycentric intersection
  // can otherwise acquire opposite roundoff signs at its edge endpoints.
  // This arithmetic bound is relative to the actual terms, never a world-
  // distance weld or an anatomical tolerance. Source-node signs are retained.
  double error=64*std::numeric_limits<double>::epsilon()*magnitude;
  return std::abs(d)<=error?0:d;
 }
 ID cut_plane(ID a,ID b,size_t plane){
  double da=plane_value(a,plane),db=plane_value(b,plane);if(da==0)return a;if(db==0)return b;if(da*db>=0)throw std::runtime_error("noncrossing scalar edge");std::set<ID> support;for(auto [i,w]:points[a].weights)support.insert(i);for(auto [i,w]:points[b].weights)support.insert(i);if(support.size()>4||support.empty())throw std::runtime_error("invalid scalar simplex");ID ta=0,tb=0;for(size_t j=0;j<planes;j++){if(plane_value(a,j)==0)ta|=ID(1)<<j;if(plane_value(b,j)==0)tb|=ID(1)<<j;}
  CutKey key;key.labels=(ta&tb)|(ID(1)<<plane);key.size=support.size();std::copy(support.begin(),support.end(),key.nodes.begin());auto found=cache.find(key);if(found!=cache.end())return found->second;std::vector<size_t> tied;for(size_t j=0;j<planes;j++)if(key.labels&(ID(1)<<j))tied.push_back(j);std::vector<std::vector<double>> A(tied.size()+1,std::vector<double>(key.size,1.));std::vector<double> rhs(tied.size()+1);rhs[0]=1;
  for(size_t row=0;row<tied.size();row++)for(size_t k=0;k<key.size;k++)A[row+1][k]=score[key.nodes[k]][tied[row]];std::vector<double>w;try{w=least_squares(A,rhs);}catch(const std::exception& error){std::cerr<<"scalar cut plane "<<plane<<"; simplex";for(size_t k=0;k<key.size;k++)std::cerr<<" "<<key.nodes[k];std::cerr<<"; tied planes";for(auto j:tied)std::cerr<<" "<<j;std::cerr<<"\n"<<std::setprecision(17);for(auto& row:A){for(double x:row)std::cerr<<x<<" ";std::cerr<<"\n";}for(ID endpoint:{a,b}){std::cerr<<"edge endpoint "<<endpoint<<"; ties "<<points[endpoint].ties<<"; plane0 "<<plane_value(endpoint,0)<<"; plane1 "<<plane_value(endpoint,1)<<"; weights";for(auto [node,weight]:points[endpoint].weights)std::cerr<<" "<<node<<":"<<weight;std::cerr<<"\n";}throw;}
  double residual=0;for(size_t row=0;row<A.size();row++){double d=-rhs[row];for(size_t k=0;k<w.size();k++)d+=A[row][k]*w[k];residual=std::max(residual,std::abs(d));}double violation=0;for(double x:w)violation=std::max(violation,std::max(-x,x-1.));max_residual=std::max(max_residual,residual);max_weight_error=std::max(max_weight_error,violation);if(residual>2e-11||violation>2e-11)throw std::runtime_error("scalar cut leaves source simplex");
  // Keep the minimal source simplex. Solving redundant incidences can leave
  // a roundoff-sized coefficient on a node whose exact weight is zero. Such
  // a node otherwise makes a later edge appear to have four independent
  // supports while it is actually on a source face. This arithmetic bound
  // is far below the independently checked 2e-11 barycentric residual gate.
  // Use the same verified source-simplex reserve as the preceding weight
  // gate. Keeping a tolerated negative coefficient as a support node is
  // inconsistent: a point on a source edge then masquerades as face-interior.
  double reserve=2e-11,total=0;
  for(double& x:w){if(std::abs(x)<=reserve){max_pruned_weight=std::max(max_pruned_weight,std::abs(x));x=0;}total+=x;}for(double& x:w)x/=total;
  for(size_t row=0;row<A.size();row++){double d=-rhs[row];for(size_t k=0;k<w.size();k++)d+=A[row][k]*w[k];if(std::abs(d)>2e-11)throw std::runtime_error("minimal-simplex replay failed");}
  Point q{};q.ties=key.labels;for(size_t k=0;k<w.size();k++)if(w[k]!=0){q.p=add(q.p,mul(points[key.nodes[k]].p,w[k]));q.weights.push_back({key.nodes[k],w[k]});}ID index=points.size();points.push_back(q);cache.emplace(key,index);return index;
 }
 ID cut(ID a,ID b,size_t owner,size_t other){
  double da=value(a,owner,other),db=value(b,owner,other);if(da==0)return a;if(db==0)return b;if(da*db>=0)throw std::runtime_error("noncrossing calculation edge");
  std::set<ID> support;for(auto [i,w]:points[a].weights)support.insert(i);for(auto [i,w]:points[b].weights)support.insert(i);
  if(support.size()>4||support.empty())throw std::runtime_error("invalid source simplex");ID ta=0,tb=0;
  for(size_t j=0;j<labels;j++){if(value(a,owner,j)==0)ta|=ID(1)<<j;if(value(b,owner,j)==0)tb|=ID(1)<<j;}
  CutKey key;key.labels=(ta&tb)|(ID(1)<<owner)|(ID(1)<<other);key.size=support.size();std::copy(support.begin(),support.end(),key.nodes.begin());auto found=cache.find(key);if(found!=cache.end())return found->second;
  std::vector<size_t> tied;for(size_t j=0;j<labels;j++)if(key.labels&(ID(1)<<j))tied.push_back(j);
  std::vector<std::vector<double>> A(tied.size(),std::vector<double>(key.size,1.));std::vector<double> rhs(tied.size());rhs[0]=1.;
  for(size_t row=1;row<tied.size();row++)for(size_t k=0;k<key.size;k++)A[row][k]=score[key.nodes[k]][tied[row]]-score[key.nodes[k]][tied[0]];
  auto w=least_squares(A,rhs);double residual=0;for(size_t row=0;row<A.size();row++){double d=-rhs[row];for(size_t k=0;k<w.size();k++)d+=A[row][k]*w[k];residual=std::max(residual,std::abs(d));}
  double violation=0;for(double x:w)violation=std::max(violation,std::max(-x,x-1.));max_residual=std::max(max_residual,residual);max_weight_error=std::max(max_weight_error,violation);
  if(residual>2e-11||violation>2e-11)throw std::runtime_error("canonical cut leaves its source simplex");
  Point q{};q.ties=key.labels;for(size_t k=0;k<w.size();k++){q.p=add(q.p,mul(points[key.nodes[k]].p,w[k]));q.weights.push_back({key.nodes[k],w[k]});}
  ID index=points.size();points.push_back(q);cache.emplace(key,index);return index;
 }
 std::vector<ID> ordered(std::set<ID> ids,V normal){
  if(ids.size()<3)throw std::runtime_error("collapsed region polygon");V center{};for(ID i:ids)center=add(center,points[i].p);center=mul(center,1./ids.size());normal=unit(normal);size_t axis=0;for(size_t k=1;k<3;k++)if(std::abs(normal[k])<std::abs(normal[axis]))axis=k;V basis{};basis[axis]=1;V u=unit(cross(normal,basis)),v=cross(normal,u);
  std::vector<std::pair<double,ID>> angles;for(ID i:ids){V d=sub(points[i].p,center);angles.push_back({std::atan2(dot(d,v),dot(d,u)),i});}std::sort(angles.begin(),angles.end());std::vector<ID> result;for(auto [a,i]:angles)result.push_back(i);auto smallest=std::min_element(result.begin(),result.end());std::rotate(result.begin(),smallest,result.end());return result;
 }
 void add_face(size_t label,F f){
  V a=points[f[0]].p,b=points[f[1]].p,c=points[f[2]].p;V n=cross(sub(b,a),sub(c,a));
  // The intersection of two closed inequalities may have only line/point
  // dimension. It is not a volumetric boundary facet. Omit exact zero area
  // only; no distance/area tolerance erases small finite anatomical tissue.
  if(dot(n,n)==0){omitted_zero_area++;return;}
  F key=f;std::sort(key.begin(),key.end());auto& bank=faces[label];auto prior=bank.find(key);
  if(prior==bank.end()){bank.emplace(key,f);return;}
  F old=prior->second;bool opposite=false;for(size_t k=0;k<3;k++)if(old[k]==f[0]&&old[(k+2)%3]==f[1]&&old[(k+1)%3]==f[2])opposite=true;
  if(!opposite)throw std::runtime_error("duplicate nonreciprocal calculation face");bank.erase(prior);
 }
 void emit_polygon(size_t label,const std::vector<ID>& polygon){
  if(polygon.size()==3){add_face(label,{polygon[0],polygon[1],polygon[2]});return;}
  auto key=polygon;std::sort(key.begin(),key.end());ID c;auto found=centroids.find(key);
  if(found!=centroids.end())c=found->second;else{Point q{};for(ID i:polygon)q.p=add(q.p,points[i].p);q.p=mul(q.p,1./polygon.size());c=points.size();points.push_back(q);centroids.emplace(key,c);}
  for(size_t i=0;i<polygon.size();i++)add_face(label,{c,polygon[i],polygon[(i+1)%polygon.size()]});
 }
 void add_polygon(size_t label,const std::vector<ID>& polygon){
  // A nested anatomical subdivision must also refine its neighbors' shared
  // facets. Split their boundary charts without adding an internal material
  // cap or changing occupied volume.
  std::vector<std::vector<ID>> polygons{polygon};
  for(ID plane:surface_refinements){std::vector<std::vector<ID>> next;
   for(const auto& face:polygons){std::vector<double> d;bool positive=false,negative=false;for(ID node:face){double value=plane_value(node,plane);d.push_back(value);positive|=value>0;negative|=value<0;}
    if(!positive||!negative){next.push_back(face);continue;}
    for(int side:{1,-1}){std::vector<ID> part;for(size_t k=0;k<face.size();k++){size_t j=(k+1)%face.size();if(side*d[k]>=0)part.push_back(face[k]);if(d[k]*d[j]<0)part.push_back(cut_plane(face[k],face[j],plane));}if(part.size()>=3)next.push_back(std::move(part));}
   }polygons=std::move(next);
  }for(const auto& face:polygons)emit_polygon(label,face);
 }
 void tetrahedron(const std::array<ID,4>& t){
  V a=sub(points[t[1]].p,points[t[0]].p),b=sub(points[t[2]].p,points[t[0]].p),c=sub(points[t[3]].p,points[t[0]].p);if(!(dot(a,cross(b,c))>0))throw std::runtime_error("nonpositive calculation tetrahedron");
  std::vector<size_t> active;for(size_t i=0;i<labels;i++){
   bool dominated=false;for(size_t j=0;j<labels&&!dominated;j++)if(i!=j){bool ge=true,gt=false;for(ID node:t){ge&=score[node][j]>=score[node][i];gt|=score[node][j]>score[node][i];}if(ge&&(gt||j<i))dominated=true;}
   if(!dominated)active.push_back(i);
  }
  for(size_t owner:active){
   std::vector<std::vector<ID>> polygons{{t[1],t[2],t[3]},{t[0],t[3],t[2]},{t[0],t[1],t[3]},{t[0],t[2],t[1]}};
   for(size_t other:active)if(owner!=other&&!polygons.empty()){
    std::set<ID> all;for(auto& f:polygons)all.insert(f.begin(),f.end());std::map<ID,double>d;bool ge=true,le=true;for(ID i:all){d[i]=value(i,owner,other);ge&=d[i]>=0;le&=d[i]<=0;}if(ge)continue;if(le){polygons.clear();break;}
    std::vector<std::vector<ID>> next;std::set<ID> cap;
    for(auto& face:polygons){std::vector<ID> polygon;for(size_t k=0;k<face.size();k++){ID a=face[k],b=face[(k+1)%face.size()];if(d[a]>=0)polygon.push_back(a);if(d[a]*d[b]<0){ID i=cut(a,b,owner,other);polygon.push_back(i);cap.insert(i);}else if(d[a]==0)cap.insert(a);}if(polygon.size()>=3)next.push_back(std::move(polygon));}
    if(cap.size()>=3){std::vector<std::vector<double>> A(3,std::vector<double>(3));std::vector<double> rhs(3);for(size_t k=0;k<3;k++){V delta=sub(points[t[k+1]].p,points[t[0]].p);for(size_t j=0;j<3;j++)A[k][j]=delta[j];rhs[k]=score[t[k+1]][owner]-score[t[k+1]][other]-score[t[0]][owner]+score[t[0]][other];}auto g=least_squares(A,rhs);next.push_back(ordered(cap,{-g[0],-g[1],-g[2]}));}polygons=std::move(next);
   }
   for(auto& polygon:polygons)add_polygon(owner,polygon);
  }
 }
 void constrained_tetrahedron(const std::array<ID,4>& t,int domain=-1){
  V a=sub(points[t[1]].p,points[t[0]].p),b=sub(points[t[2]].p,points[t[0]].p),c=sub(points[t[3]].p,points[t[0]].p);if(!(dot(a,cross(b,c))>0))throw std::runtime_error("nonpositive constrained tetrahedron");
  for(size_t owner=0;owner<labels;owner++){
   if(domain>=0&&cell_domain[owner]!=ID(domain))continue;
   bool empty=false;for(auto [plane,side]:constraints[owner]){bool allzero=true,le=true,lt=false;for(ID i:t){double d=side*score[i][plane];allzero&=d==0;le&=d<=0;lt|=d<0;}if((allzero&&side!=zero_side[plane])||(le&&lt)){empty=true;break;}}if(empty)continue;
   std::vector<std::vector<ID>> polygons{{t[1],t[2],t[3]},{t[0],t[3],t[2]},{t[0],t[1],t[3]},{t[0],t[2],t[1]}};
   for(auto [plane,side]:constraints[owner]){
    if(polygons.empty())break;std::set<ID> all;for(auto& f:polygons)all.insert(f.begin(),f.end());std::map<ID,double>d;bool ge=true,le=true;for(ID i:all){d[i]=side*plane_value(i,plane);ge&=d[i]>=0;le&=d[i]<=0;}if(ge)continue;if(le){polygons.clear();break;}
    std::vector<std::vector<ID>> next;std::set<ID> cap;
    for(auto& face:polygons){std::vector<ID> polygon;for(size_t k=0;k<face.size();k++){ID a=face[k],b=face[(k+1)%face.size()];if(d[a]>=0)polygon.push_back(a);if(d[a]*d[b]<0){ID i=cut_plane(a,b,plane);polygon.push_back(i);cap.insert(i);}else if(d[a]==0)cap.insert(a);}if(polygon.size()>=3)next.push_back(std::move(polygon));}
    if(cap.size()>=3){V a=sub(points[t[1]].p,points[t[0]].p),b=sub(points[t[2]].p,points[t[0]].p),c=sub(points[t[3]].p,points[t[0]].p);double r0=score[t[1]][plane]-score[t[0]][plane],r1=score[t[2]][plane]-score[t[0]][plane],r2=score[t[3]][plane]-score[t[0]][plane];V g=add(add(mul(cross(b,c),r0),mul(cross(c,a),r1)),mul(cross(a,b),r2));next.push_back(ordered(cap,mul(g,-double(side))));}polygons=std::move(next);
   }
   for(auto& polygon:polygons)add_polygon(owner,polygon);
  }
 }
};
int main(int argc,char** argv){
 auto start=std::chrono::steady_clock::now();try{
  if(argc!=4)throw std::runtime_error("input output finite-seconds required");double seconds=std::stod(argv[3]);if(!(seconds>0&&seconds<=300))throw std::runtime_error("invalid deadline");std::ifstream in(argv[1],std::ios::binary);char magic[8];in.read(magic,8);std::string mode(magic,8);bool material_domain=mode=="ATLPDC01"||mode=="ATLPDR01";bool refine_surfaces=mode=="ATLPLR01"||mode=="ATLPDR01";bool constrained=material_domain||refine_surfaces||mode=="ATLPLC01";if(!constrained&&mode!="ATLTET01")throw std::runtime_error("invalid input header");ID np=read<ID>(in),nt=read<ID>(in),cols=read<ID>(in),nc=constrained?read<ID>(in):cols;if(np>3000000||nt>2000000||nc<2||nc>32||cols>64)throw std::runtime_error("explicit calculation budget exceeded");
  Engine e;e.labels=nc;e.planes=constrained?cols:0;e.faces.resize(nc);e.points.reserve(np);e.score.resize(np,std::vector<double>(cols));for(ID i=0;i<np;i++){Point q{};for(double& x:q.p){x=read<double>(in);if(!std::isfinite(x))throw std::runtime_error("nonfinite position");}q.weights.push_back({i,1.});e.points.push_back(q);}std::vector<std::array<ID,4>> tets(nt);for(auto& t:tets)for(ID& i:t){i=read<ID>(in);if(i>=np)throw std::runtime_error("invalid source index");}for(auto& row:e.score)for(double& x:row){x=read<double>(in);if(!std::isfinite(x))throw std::runtime_error("nonfinite score");}
  std::vector<ID> tet_domain;
  if(constrained){for(ID i=0;i<cols;i++){int64_t side=read<int64_t>(in);if(side!=1&&side!=-1)throw std::runtime_error("invalid zero-side priority");e.zero_side.push_back(side);}if(material_domain){tet_domain.resize(nt);for(ID& d:tet_domain){d=read<ID>(in);if(d>63)throw std::runtime_error("invalid tetrahedral material domain");}e.cell_domain.resize(nc);for(ID& d:e.cell_domain){d=read<ID>(in);if(d>63)throw std::runtime_error("invalid cell material domain");}}if(refine_surfaces){ID count=read<ID>(in);if(count>cols)throw std::runtime_error("surface refinement budget");for(ID i=0;i<count;i++){ID plane=read<ID>(in);if(plane>=cols)throw std::runtime_error("invalid refinement plane");e.surface_refinements.push_back(plane);}}e.constraints.resize(nc);for(auto& cell:e.constraints){ID refs=read<ID>(in);if(refs>cols)throw std::runtime_error("constraint budget");for(ID j=0;j<refs;j++){ID plane=read<ID>(in);int64_t side=read<int64_t>(in);if(plane>=cols||(side!=1&&side!=-1))throw std::runtime_error("invalid cell constraint");cell.push_back({plane,side});}}}
  for(ID i=0;i<nt;i++){if(i%256==0){double elapsed=std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count();if(elapsed>seconds)throw std::runtime_error("explicit native partition deadline");if(i&&i%50000==0)std::cerr<<"partition "<<i<<" / "<<nt<<"; "<<elapsed<<" s\n";}if(constrained)e.constrained_tetrahedron(tets[i],material_domain?int(tet_domain[i]):-1);else e.tetrahedron(tets[i]);}
  for(auto& cell:e.faces)for(auto& [key,f]:cell){V a=e.points[f[0]].p,b=e.points[f[1]].p,c=e.points[f[2]].p;if(!(dot(cross(sub(b,a),sub(c,a)),cross(sub(b,a),sub(c,a)))>0))throw std::runtime_error("degenerate output facet");}
  std::string tmp=std::string(argv[2])+".tmp";std::ofstream out(tmp,std::ios::binary);out.write("ATLPAR01",8);write(out,ID(e.points.size()));write(out,nc);for(auto& q:e.points)for(double x:q.p)write(out,x);for(auto& cell:e.faces){write(out,ID(cell.size()));for(auto& [key,f]:cell)for(ID i:f)write(out,i);}out.close();if(!out)throw std::runtime_error("partition output failed");if(std::rename(tmp.c_str(),argv[2])!=0)throw std::runtime_error("atomic output rename failed");double elapsed=std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count();std::cout<<"{\"passed\":true,\"tetrahedra\":"<<nt<<",\"positions\":"<<e.points.size()<<",\"labels\":"<<nc<<",\"maximum_score_residual\":"<<e.max_residual<<",\"maximum_weight_violation\":"<<e.max_weight_error<<",\"elapsed_seconds\":"<<elapsed<<"}\n";return 0;
 }catch(const std::exception& error){std::cerr<<error.what()<<"\n";return 1;}
}
