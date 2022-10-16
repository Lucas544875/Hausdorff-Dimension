let fragmentShader =`
precision mediump float;
uniform float time;
uniform vec2  resolution;
uniform vec3  cDir;
uniform vec3  cPos;

const float PI = 3.14159265;
const float PI2 = PI*2.0;
const float E = 2.71828182;
const float INFINITY = 1.e20;
const float FOV = 30.0 * 0.5 * PI / 180.0;//field of view
const vec3 LightDir = normalize(vec3(2.0,1.0,1.0));
const int Iteration =128;
const int MAX_REFRECT = 2;

struct rayobj{
  vec3  rPos;     //レイの場所
  vec3  direction;//方向
  float distance; //距離関数の返り値
  float len;      //出発点からの距離
  float iterate;  //レイマーチの反復回数
  int   objectID;  //オブジェクトID
  int   material; //マテリアルID
  vec3  normal;   //法線ベクトル
  vec3  fragColor;//色
};

struct effectConfig{
  bool reflect;    //反射
  bool ambient;    //アンビエント
  bool specular;   //ハイライト(鏡面反射)
  bool diffuse;    //拡散光
  bool incandescence;//白熱光
  bool shadow;     //ソフトシャドウ
  bool globallight;//大域照明
  bool grow;       //グロー
  bool fog;        //霧
  bool gamma;      //ガンマ補正
};

const effectConfig effect = effectConfig(
  false, //反射
  true,  //アンビエント
  false, //ハイライト(鏡面反射)
  true, //拡散光
  false,  //白熱光
  true,  //ソフトシャドウ
  false, //大域照明
  false, //グロー
  false,  //霧
  true   //ガンマ補正
);

struct dfstruct{
  float dist;
  int   id;
};


//quaternion
vec4 times(vec4 q1,vec4 q2){
  return vec4 (
    q1[0]*q2[0] - q1[1]*q2[1] - q1[2]*q2[2] - q1[3]*q2[3],
    q1[0]*q2[1] + q1[1]*q2[0] + q1[2]*q2[3] - q1[3]*q2[2],
    q1[0]*q2[2] - q1[1]*q2[3] + q1[2]*q2[0] + q1[3]*q2[1],
    q1[0]*q2[3] + q1[1]*q2[2] - q1[2]*q2[1] + q1[3]*q2[0]
  );
}

vec4 inverse(vec4 q){
  return vec4(q[0],-q[1],-q[2],-q[3]);
}

vec4 rotation(float theta,vec3 v){
  float c = cos(theta/2.0);
  float s = sin(theta/2.0);
  return normalize(vec4(c,v.x*s,v.y*s,v.z*s));
}

vec4 turn(vec4 v,vec4 rot){
  return times(times(rot,v),inverse(rot));
}

vec3 turn(vec3 v,vec3 pole,float theta){//vをpole中心に角theta回転
  vec4 rot = rotation(theta,pole);
  return turn(vec4(0,v),rot).yzw;
}

float angle(vec3 v, vec3 w){
  vec3 nv = normalize(v);
  vec3 nw = normalize(w);
  return acos(dot(nv,nw));
}

vec3 hsv(float h, float s, float v) {
  // h: 0.0 - 2PI, s: 0.0 - 1.0, v: 0.0 - 1.0, 円柱モデル
  return ((clamp(abs(fract(mod(h,2.0*PI)+vec3(0,2,1)/3.)*6.-3.)-1.,0.,1.)-1.)*s+1.)*v;
}

float manhattan (vec3 p,vec3 q){
  return abs(p.x-q.x)+abs(p.y-q.y)+abs(p.z-q.z);
}

float chebyshev (vec3 p,vec3 q){
  return max(max(abs(p.x-q.x),abs(p.y-q.y)),abs(p.z-q.z));
}

vec3 Hadamard(vec3 v,vec3 w){ //アダマール積
  return vec3(
    v.x * w.x,
    v.y * w.y,
    v.z * w.z
  );
}

vec3 gammainv(vec3 p){
  return pow(p,vec3(1.0/2.2));
}

mat2 rot(float a) {
  float c = cos(a), s = sin(a);
  return mat2(c,s,-s,c);
}

mat3 transpose(mat3 m) {
  return mat3(m[0][0], m[1][0], m[2][0],
              m[0][1], m[1][1], m[2][1],
              m[0][2], m[1][2], m[2][2]);
}

float det(mat2 matrix) {
  return matrix[0].x * matrix[1].y - matrix[0].y * matrix[1].x;
}

mat3 inverse(mat3 matrix) {
  vec3 row0 = matrix[0];
  vec3 row1 = matrix[1];
  vec3 row2 = matrix[2];

  vec3 minors0 = vec3(
    det(mat2(row1.y, row1.z, row2.y, row2.z)),
    det(mat2(row1.z, row1.x, row2.z, row2.x)),
    det(mat2(row1.x, row1.y, row2.x, row2.y))
  );
  vec3 minors1 = vec3(
    det(mat2(row2.y, row2.z, row0.y, row0.z)),
    det(mat2(row2.z, row2.x, row0.z, row0.x)),
    det(mat2(row2.x, row2.y, row0.x, row0.y))
  );
  vec3 minors2 = vec3(
    det(mat2(row0.y, row0.z, row1.y, row1.z)),
    det(mat2(row0.z, row0.x, row1.z, row1.x)),
    det(mat2(row0.x, row0.y, row1.x, row1.y))
  );

  mat3 adj = transpose(mat3(minors0, minors1, minors2));

  return (1.0 / dot(row0, minors0)) * adj;
}
vec3 mix3 (vec3 v1, vec3 v2, vec3 v3, float k){
  float c1 = max(1.0-2.0*k,0.0);
  float c2 = 1.0-2.0*abs(k-0.5);
  float c3 = max(2.0*k-1.0,0.0);
  return c1*v1 + c2*v2 + c3*v3;
}


dfstruct dfmeta(dfstruct df1, dfstruct df2,float k){ //メタボール風の結合
  float distmin, distmax;
  int id;
  if (df1.dist < df2.dist){
    distmin = df1.dist;
    distmax = df2.dist;
    id = df1.id;
  }else{
    distmin = df2.dist;
    distmax = df1.dist;
    id = df2.id;
  }
  float h = 1.0 + exp(-k *(distmax-distmin));
  return dfstruct(distmin -log(h) / k, id);
}

dfstruct dfmax(dfstruct df1, dfstruct df2){ //共通部分
  if (df1.dist < df2.dist){
    return df2;
  }else{
    return df1;
  }
}

dfstruct dfmin(dfstruct df1, dfstruct df2){//和集合
  if (df1.dist < df2.dist){
    return df1;
  }else{
    return df2;
  }
}

vec2 pmod2d(vec2 p, float r,float section) {
  float a = atan(p.x, p.y) + PI/r+section;
  float n = PI2 / r;
  a = floor(a/n)*n ;
  return p*rot(-a);
}

vec3 pmod(vec3 z, vec3 center, vec3 direction, int n, float section){
  vec3 cz = z - center;
  vec3 pole = cross(vec3(0,0,1),direction);
  float theta = angle(vec3(0,0,1),direction);
  vec3 tz = turn(cz,pole,-theta);
  vec3 zz = vec3(pmod2d(tz.xy,float(n),section),tz.z);
  return turn(zz,pole,theta) + center;
}

vec3 wipe(vec3 z, vec3 center, vec3 direction, int n, float section){//結局よくわからん
  vec3 cz = z - center;
  vec3 axes = cross(direction,vec3(0,0,1));
  float theta = angle(axes,cz) + section;
  float shift = floor(theta*float(n)/PI2)*PI2/float(n);
  return turn(cz,direction,shift)+center;
}

//primitives
float sphere(vec3 z,vec3 center,float radius){
  return length(z-center)-radius;
}


float plane(vec3 z,vec3 normal,float offset){
	return dot(z,normalize(normal)) - offset;
}

float plane1(vec3 z){//plane
  return plane(z,normalize(vec3(0.0,0.0,1.0)),0.5);
}

void sphereFold(inout vec3 z, inout float dz, float minRadius2, float fixedRadius2) {
	float r2 = dot(z,z);
	if (r2<minRadius2) { 
		// linear inner scaling
		float temp = fixedRadius2/(minRadius2);
		z *= temp;
		dz*= temp;
	} else if (r2<fixedRadius2) { 
		// this is the actual sphere inversion
		float temp =fixedRadius2/r2;
		z *= temp;
		dz*= temp;
	}
}

void boxFold(inout vec3 z, inout float dz, float foldingLimit) {
	z = clamp(z, -foldingLimit, foldingLimit) * 2.0 - z;
}

float mandelBox(vec3 z, float Scale, float foldingLimit, float minRadius2, float fixedRadius2){
	vec3 offset = z;
	float dr = 1.0;
	for (int n = 0; n < 16; n++) {
		boxFold(z,dr,foldingLimit);       // Reflect
		sphereFold(z,dr,minRadius2,fixedRadius2);    // Sphere Inversion
    z=Scale*z + offset;  // Scale & Translate
    dr = dr*abs(Scale)+1.0;
	}
	float r = length(z);
	return r/abs(dr);
}
float tofu(vec3 z){
  float Scale = 1.9 ;//定数
  float foldingLimit=0.6;//定数
	float minRadius2=0.60;//定数
  float fixedRadius2=2.65;//定数
	return mandelBox(z, Scale, foldingLimit, minRadius2, fixedRadius2);
}
float kado(vec3 z){
  float Scale = -2.18 ;//定数
  float foldingLimit=1.14;//定数
	float minRadius2=0.60;//定数
  float fixedRadius2=2.65;//定数
	return mandelBox(z, Scale, foldingLimit, minRadius2, fixedRadius2);
}

float sdCross(vec3 p, float c) {
	p = abs(p);
	float dxy = max(p.x, p.y);
	float dyz = max(p.y, p.z);
	float dxz = max(p.x, p.z);
	return min(dxy, min(dyz, dxz)) - c;
}

float sdBox(vec3 p, vec3 b) {
	p = abs(p) - b;
	return length(max(p, 0.0)) + min(max(p.x, max(p.y, p.z)), 0.0);
}

float _mengerSponge(vec3 p, float scale, float width) {
	float d = sdBox(p, vec3(1.0));
	float s = 1.0;
	for (int i = 0; i < 8; i++) {
		vec3 a = mod(p * s, 2.0) - 1.0;
		s *= scale;
		vec3 r = 1.0 - scale * abs(a);
		float c = sdCross(r, width) / s;
		d = max(d, c);
	}
	return d;
}

float mengerSponge(vec3 p) {
	float scale = 3.0;
	float width = 1.0;
	return _mengerSponge(p,scale,width);
}

float pseudoKleinian(vec3 p) {
	vec3 csize = vec3(0.90756, 0.92436, 0.90756);
	float size = 1.2;
	vec3 c = vec3(0.0);
	float defactor = 1.0;
	vec3 ap = p + 1.0;
	for (int i = 0; i < 10; i++) {
		ap = p;
		p = 2.0 * clamp(p, -csize, csize) - p;
		float r2 = dot(p, p);
		float k = max(size / r2, 1.0);
		p *= k;
		defactor *= k;
		p += c;
	}
	float r = abs(0.5 * p.z / defactor);
	return r;
}
float shiftKeinian(vec3 p){
	vec3 o = vec3(11.0,-1.0,-1.0);
	return pseudoKleinian(p - o);
}

float gasket(vec3 z){
	const float isq3 = inversesqrt(3.0);
	vec3 a1 = vec3(0,0,3.0/2.0);
	vec3 a2 = vec3(0,1,0);
	vec3 a3 = vec3(2.0/3.0*isq3,-1.0/3.0*isq3,0);
	vec3 a4 = vec3(-2.0/3.0*isq3,-1.0/3.0*isq3,0);
	vec3 c;
	float dist, d;
  float Scale=2.0;
  const int ite=50;
	for (int i=0;i < ite;i++) {
    c = a1; dist = length(z-a1);
    d = length(z-a2); if (d < dist) { c = a2; dist=d; }
    d = length(z-a3); if (d < dist) { c = a3; dist=d; }
    d = length(z-a4); if (d < dist) { c = a4; dist=d; }
    z = Scale*z-c*(Scale-1.0);
	}
	return length(z) * pow(Scale, float(-ite));
}

float lenPtoL(vec3 p,vec3 l, vec3 dir){
	vec3 ndir = normalize(dir);
	return length(l-p - dot(ndir,l-p)*ndir);
}
float triangle(vec3 z, vec3 p1, vec3 p2, vec3 p3){
	vec3 p1z = z-p1;
	vec3 p1p2 = p2 - p1;
	vec3 p1p3 = p3 - p1;
	vec3 p2p3 = p3 - p2;
	vec3 normal = normalize(cross(p1p2,p1p3));
	mat3 cmat = mat3(p1p2,p1p3,normal);
	vec3 c = inverse(cmat)*p1z;
	if(c.x > 0.0 && c.y > 0.0 && c.x+c.y < 1.0){
		return abs(dot(normal,p1z));
	}else if(c.y<0.0 && c.x>0.0 && c.x+c.y<1.0){
		return lenPtoL(z,p1,p1p2);
	}else if(c.x<0.0 && c.y>0.0 && c.y+c.x<1.0){
		return lenPtoL(z,p1,p1p3);
	}else if(c.x+c.y>1.0 && abs(c.x-c.y)<1.0){
		return lenPtoL(z,p2,p2p3);
	}else{
		float d1 = length(z-p1);
		float d2 = length(z-p2);
		float d3 = length(z-p3);
		return min(min(d1,d2),d3);
	}
}

float octahedron(vec3 z){
  vec3 pz = pmod(z,vec3(0),vec3(1,0,0),2,PI/2.0);
  vec3 ppz = pmod(pz,vec3(0),vec3(0,0,1),4,-PI/4.0);
  vec3 p1 = vec3(1,0,0);
  vec3 p2 = vec3(0,1,0);
  vec3 p3 = vec3(0,0,1);
  return triangle(ppz,p1,p2,p3)-0.01;
}


float box(vec3 p){
  p.xy*=rot(3.0/4.0*PI -time*0.05);
  p.xz*=rot(PI/2.0-atan(sqrt(3.0),sqrt(2.0)));
  p.yz*=rot(PI/4.0);
  p=abs(p);
  p-=1.0;
  if(p.x<p.y)p.xy=p.yx;
  if(p.x<p.z)p.xz=p.zx;
  if(p.y<p.z)p.yz=p.zy;
  return max(abs(p.x),abs(p.y)) -0.12;
}
float floor1(vec3 z){
  return plane(z,vec3(0,0,1),-1.95);
}
dfstruct distanceFunction(vec3 z){
  z=z +vec3(-2,0,0);
  dfstruct box = dfstruct(box(z),0);
  dfstruct plane = dfstruct(floor1(z),1);
  return dfmin(box,plane);
}
dfstruct depthFunction(vec3 z){
  z=z +vec3(-2,0,0);
  dfstruct box = dfstruct(box(z),0);
  return box;
}

//マテリアルの設定
const int SAIHATE = 0;
const int CYAN = 1;
const int WHITE = 2;
const int GRID = 3;
const int MANDEL = 4;
const int BROWN = 5;
const int NORMAL = 6;
const int METAL = 7;
const int KADO = 8;
const int MAT = 9;
const int LESSSTEP = 97;
const int DEBUG = 98;
const int ERROR = 99;

//マテリアルの設定
int materialOf(int objectID){
  if (objectID == 0){
    return KADO;
  }else if (objectID == 1){
    return WHITE;
  }else if (objectID == 2){
    return WHITE;
  }else if (objectID == 98){
    return SAIHATE;
  }else if (objectID == 99){
    return LESSSTEP;
  }else{
    return ERROR;
  }
}

vec3 normal(vec3 p){
  float d = 0.005;
  return normalize(vec3(
    distanceFunction(p + vec3(  d, 0.0, 0.0)).dist - distanceFunction(p + vec3( -d, 0.0, 0.0)).dist,
    distanceFunction(p + vec3(0.0,   d, 0.0)).dist - distanceFunction(p + vec3(0.0,  -d, 0.0)).dist,
    distanceFunction(p + vec3(0.0, 0.0,   d)).dist - distanceFunction(p + vec3(0.0, 0.0,  -d)).dist
  ));
}


vec3 gridCol(vec3 rPos){
  return mix(vec3(0.3),vec3(step(fract(2.0*rPos.x),0.05),step(fract(2.0*rPos.y),0.05),step(fract(2.0*rPos.z),0.05)),0.5);
}

vec3 debugCol(vec3 rPos){
  return fract(rPos);
}

vec3 kadoCol(vec3 rPos){
  float k = dot (rPos,vec3(0,1,1))/8.0 +0.5;
  vec3 col1 = gammainv(vec3(0.007, 0.313, 0.772));
  vec3 col2 = gammainv(vec3(0.831, 0.247, 0.552));
  return mix(col1,col2,k);
}

vec3 normalCol(vec3 rPos){
  return abs(normal(rPos));
}

vec3 color(rayobj ray){
  if (ray.material == GRID){
    return gridCol(ray.rPos);
  }else if (ray.material == WHITE){
    return vec3(1.0,1.0,1.0);
  }else if (ray.material == DEBUG){
    return debugCol(ray.rPos);
  }else if (ray.material == MANDEL){
    return kadoCol(ray.rPos);
  }else if (ray.material == BROWN){
    return vec3(0.454, 0.301, 0.211);
  }else if (ray.material == NORMAL){
    return normalCol(ray.rPos);
  }else if (ray.material == METAL){
    return vec3(0.7);
  }else if (ray.material == KADO){
    return kadoCol(ray.rPos);
  }else if (ray.material == MAT){
    return vec3(0.960,0.95,0.92);
  }else if (ray.material == SAIHATE || ray.material == LESSSTEP){
    float k = max(0.0,dot(normalize(ray.direction),vec3(0,0,1)));
    vec3 c1 = vec3(1.0);
    vec3 c2 = vec3(0.584, 0.752, 0.925);
    return mix(c1,c2,smoothstep(0.0,0.5,k));
    //return vec3(160.0,216.0,239.0)/256.0;
  }else{
    return vec3(1.0,0.0,0.0);
  }
}

float refrectance(int material){
  if (material == CYAN){
    return 0.1;
  }else if (material == WHITE){
    return 0.6;
  }else if (material == DEBUG){
    return 0.3;
  }else if (material == GRID){
    return 0.3;
  }else if (material == MANDEL){
    return 0.3;
  }else if (material == NORMAL){
    return 0.4;
  }else if (material == METAL){
    return 0.6;
  }else{
    return 0.0;
  }
}


void raymarch(inout rayobj ray){
  for(int i = 0; i < Iteration; i++){
    dfstruct df = distanceFunction(ray.rPos);
    ray.distance = df.dist;
    if(ray.distance < 0.001){
      ray.normal = normal(ray.rPos);
      ray.objectID = df.id;
      ray.iterate = float(i)/float(Iteration);
      return;
    }
    ray.len += ray.distance;
    if(ray.len > 100.0){
      ray.objectID = 98;
      ray.iterate = float(i)/float(Iteration);
      return;
    }
    ray.rPos += ray.distance * ray.direction;
  }
  ray.objectID = 99;
  ray.iterate = 1.0;
}

void trick(inout rayobj ray){
  rayobj probe = ray;
  const float overstep = 2.0;
  probe.rPos += overstep * probe.direction;
  probe.len += overstep;

  for(int i = 0; i < Iteration; i++){
    dfstruct df = depthFunction(probe.rPos);
    probe.distance = df.dist;
    if(probe.distance<0.0){
      probe.rPos += overstep * probe.direction;
      probe.len += overstep;
      continue;
    }else if(probe.distance < 0.001){
      probe.normal = normal(probe.rPos);
      probe.objectID = df.id;
      probe.iterate = float(i)/float(Iteration);
      ray = probe;
      return;
    }
    probe.len += probe.distance;
    if(probe.len > 100.0){
      return;
    }
    probe.rPos += probe.distance * probe.direction;
  }
}

//ライティング
void ambientFunc(inout rayobj ray){//アンビエント
  vec3 baseColor = color(ray);
  vec3 ambColor = vec3(1.0);
  float ambIntensity =  0.6;
  ray.fragColor += ambIntensity * Hadamard(baseColor,ambColor);
  ray.fragColor = clamp(ray.fragColor,0.0,1.0);
}

void specularFunc(inout rayobj ray){//鏡面反射
  float t = -dot(ray.direction,ray.normal);
  vec3 reflection=ray.direction+2.0*t*ray.normal;
  float x = dot(reflection,LightDir);
  float specular=1.0/(50.0*(1.001-clamp(x,0.0,1.0)));
  ray.fragColor = clamp(ray.fragColor+specular,0.0,1.0);
}

void diffuseFunc(inout rayobj ray){//拡散光
  vec3 color = color(ray);
  vec3 lightColor = vec3(1.0);//(0.741, 0.741, 0.717);
  float diffIntensity = 1.1;
  float diffuse = max(0.0,dot(LightDir, ray.normal));
  ray.fragColor += diffIntensity * diffuse * Hadamard(color,lightColor);
  ray.fragColor = clamp(ray.fragColor,0.0,1.0);
}

void _incandescenceFunc(inout rayobj ray, vec3 incandescenceColor, vec3 incCenter, float incRadius, float incIntensity){ 
  vec3 color = pow(max((1.0 - (length(incCenter - ray.rPos)/incRadius)),0.0),2.0) * incIntensity * incandescenceColor;
  ray.fragColor += color;
  ray.fragColor = clamp(ray.fragColor,0.0,1.0);
}

void incandescenceFunc(inout rayobj ray){ //白熱光
  vec3 incandescenceColor = vec3(0.000, 0.447, 0.737);
  vec3 incCenter = vec3(0.0);
  float incRadius = 13.0;
  float incIntensity = 1.0;
  _incandescenceFunc(ray, incandescenceColor, incCenter, incRadius, incIntensity);

}

const float shadowCoef = 0.7;
void shadowFunc(inout rayobj ray){
  if (dot(ray.normal, LightDir)<0.0){return;}
  float h = 0.0;
  float c = 0.0;
  float r = 1.0;
  for(float t = 0.0; t < 50.0; t++){
    h = distanceFunction(ray.rPos + ray.normal*0.001 + LightDir * c).dist;
    if(h < 0.001){
      ray.fragColor *= shadowCoef;
      return;
    }
    r = min(r, h * 200.0 / c);
    c += h;
  }
  ray.fragColor *= mix(shadowCoef, 1.0, r);
  return;
}

void globallightFunc(inout rayobj ray){//大域照明
  vec3 origin = ray.rPos+ray.normal*0.001;
  rayobj ray2 = rayobj(origin,ray.normal,0.0,0.0,0.0,99,0,vec3(0.0),vec3(0.0));
  raymarch(ray2);
  float near = 0.10;
  ray.fragColor *= clamp(min(near,ray2.len)/near,0.0,1.0);
}

void skysphereFunc(inout rayobj ray){//天球
  if (ray.objectID == 98){
    ray.fragColor += color(ray);
  }
}

void lessStepFunc(inout rayobj ray){
  if (ray.objectID == 99){
    ray.fragColor += color(ray);
  }
}

const float growIntencity = 1.0;
void growFunc(inout rayobj ray){//グロー
  float coef = smoothstep(0.0,1.0,ray.iterate);
  const vec3 growCol = vec3(1.0);
  vec3 grow = growIntencity * coef * growCol;
  ray.fragColor += grow;
}

const float minRadius = 10.0;
const float maxRadius = 30.0;
void fogFunc(inout rayobj ray){//霧
  rayobj ray2 = ray;
  ray2.material = SAIHATE;
  vec3 fogColor = color(ray2);
  float fogcoef = clamp((ray.len-minRadius)/(maxRadius-minRadius),0.0,1.0);
  ray.fragColor = mix(ray.fragColor, fogColor, fogcoef);
}

void gammaFunc(inout rayobj ray){//ガンマ補正
  ray.fragColor=pow(ray.fragColor,vec3(2.2));
}

void reflectFunc(inout rayobj ray){//反射
  rayobj rays[MAX_REFRECT+1];
  rays[0] = ray;
  int escape = MAX_REFRECT;
  for (int i = 0;i<MAX_REFRECT;i++){
    float dot = -dot(rays[i].direction,rays[i].normal);
    vec3 direction=rays[i].direction+2.0*dot*rays[i].normal;//refrect
    rays[i+1] = rayobj(rays[i].rPos+rays[i].normal*0.001,direction,0.0,0.0,0.0,99,0,vec3(0.0),vec3(0.0));
    raymarch(rays[i+1]);
    rays[i+1].material = materialOf(rays[i+1].objectID);

    if(abs(rays[i].distance) >= 0.001){//脱出
      escape = i;
      break;
    }
  }

  for (int i = MAX_REFRECT;i >= 1;i--){
    if (i>escape){continue;}

    if(abs(ray.distance) < 0.001){//物体表面にいる場合
      if(effect.ambient){
        ambientFunc(rays[i]);
      }
      if (effect.specular){
        specularFunc(rays[i]);
      }
      if (effect.diffuse){
        diffuseFunc(rays[i]);
      }
      if (effect.incandescence){
        incandescenceFunc(rays[i]);
      }
      if (effect.shadow){
        shadowFunc(rays[i]);
      }
      if (effect.globallight){
        globallightFunc(rays[i]);
      }
    }else{//描写範囲外 or ステップ数不足
      skysphereFunc(rays[i]);
    }
    //全体
    if (effect.grow){
      growFunc(rays[i]);
    }
    if (effect.fog){
      fogFunc(rays[i]);
    }

    float refrectance = refrectance(rays[i-1].material);
    rays[i-1].fragColor += refrectance*rays[i].fragColor;
  }
  ray.fragColor += rays[0].fragColor;
}

void main(void){
  // fragment position
  vec2 p = (gl_FragCoord.xy * 2.0 - resolution) / min(resolution.x, resolution.y);
  //fix:真上真下が見えない
  vec3 xAxes = normalize(cross(cDir,vec3(0.0,0.0,1.0)));
  vec3 yAxes = normalize(-cross(cDir,xAxes));
  vec4 rot = normalize(times(rotation(-p.x*FOV,yAxes),rotation(p.y*FOV,xAxes)));
  vec3 direction = normalize(turn(vec4(0,cDir),rot).yzw);

  //レイの定義と移動
  rayobj ray = rayobj(cPos,direction,0.0,0.0,0.0,99,0,vec3(0.0),vec3(0.0));
  raymarch(ray);
  trick(ray);//錯視
  ray.material = materialOf(ray.objectID);

  //エフェクト
  if(abs(ray.distance) < 0.001){//物体表面にいる場合
    if (effect.reflect){
      reflectFunc(ray);
    }
    if(effect.ambient){
      ambientFunc(ray);
    }
    if (effect.specular){
      specularFunc(ray);
    }
    if (effect.diffuse){
      diffuseFunc(ray);
    }
    if (effect.incandescence){
      incandescenceFunc(ray);
    }
    if (effect.shadow){
      shadowFunc(ray);
    }
    if (effect.globallight){
      globallightFunc(ray);
    }
  }else{//描写範囲外 or ステップ数不足
    skysphereFunc(ray);
    lessStepFunc(ray);
  }
  //全体
  if (effect.grow){
    growFunc(ray);
  }
  if (effect.fog){
    fogFunc(ray);
  }
  if (effect.gamma){
    gammaFunc(ray);
  }
  gl_FragColor = vec4(ray.fragColor,1.0);
}
`