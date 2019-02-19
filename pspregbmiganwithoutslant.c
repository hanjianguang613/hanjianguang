/***********************************************************************************************/ 
/* pspregbmiganwithoutslant: Prestack Gaussian-beam depth migration in anisotropic media(PS-wave) without local slant stack*/
/* Note: The program developed from SU program.    */
/***********************************************************************************************/ 

#include "gb.h"
#include "segy.h"

/* Ray types */
/* one step along ray */
typedef struct RayStepStruct {
	float t;		/* time */
	float x,z;		/* x,z coordinates */
	float q1,p1,q2,p2;	/* Cerveny's dynamic ray tracing solution */
	int kmah;		/* KMAH index */
	float c,s;		/* cos(angle) and sin(angle) */
	float v,dvdx,dvdz;	/* velocity and its derivatives */
} RayStep;

/* one ray */
typedef struct RayStruct {
	int nrs;		/* number of ray steps */
	RayStep *rs;		/* array[nrs] of ray steps */
	int nc;			/* number of circles */
	int ic;			/* index of circle containing nearest step */
	void *c;		/* array[nc] of circles */
} Ray;

/*---------------------------------------------------------------------*/
/* filtered complex beam data as a function of real and imaginary time */
typedef struct BeamDataStruct {
	int ntr;		/* number of real time samples */
	float dtr;		/* real time sampling interval */
	float ftr;		/* first real time sample */
	int nti;		/* number of imaginary time samples */
	float dti;		/* imaginary time sampling interval */
	float fti;		/* first imaginary time sample */
	complex **cf;		/* array[nti][ntr] of complex data */
} BeamData;

/* one cell in which to linearly interpolate complex time and amplitude */
typedef struct CellStruct {
	int live;	/* random number used to denote a live cell */
	int dead;	/* random number used to denote a dead cell */
	float tr;	/* real part of traveltime */
	float ti;	/* imaginary part of traveltime */
	float ar;	/* real part of amplitude */
	float ai;	/* imaginary part of amplitude */
} Cell;

/* structure containing information used to set and fill cells for beams*/
typedef struct CellsStruct {
	int nt;			/* number of time samples */
	float dt;		/* time sampling interval */
	float ft;		/* first time sample */
	int lx;			/* number of x samples per cell */
	int mx;			/* number of x cells */
	int nx;			/* number of x samples in velocity model */
	int snx;        /* number of x samples in single shot record */
	float dx;		/* x sampling interval */
	float fx;		/* first x sample */
	int lz;			/* number of z samples per cell */
	int mz;			/* number of z cells */
	int nz;			/* number of z samples */
	float dz;		/* z sampling interval */
	float fz;		/* first z sample */
	int live;		/* random number used to denote a live cell */
	int dead;		/* random number used to denote a dead cell */
	float wmin;		/* minimum (reference) frequency */
	float lmin;		/* minimum beamwidth for frequency wmin */
	Cell **cell;	/* cell array[mx][mz] */
	Ray *ray;		/* ray */
	BeamData *bd;	/* complex beam data as a function of complex time */
	float **g;		/* array[nx][nz] containing g(x,z) */
} Cells;

/* structure containing information used to set and fill cells for shot ray */
typedef struct ShotCellsStruct {
	int nt;			/* number of time samples */
	float dt;		/* time sampling interval */
	float ft;		/* first time sample */
	int lx;			/* number of x samples per cell */
	int mx;			/* number of x cells */
	int nx;			/* number of x samples in single shot record */
	float dx;		/* x sampling interval */
	float fx;		/* first x sample */
	int lz;			/* number of z samples per cell */
	int mz;			/* number of z cells */
	int nz;			/* number of z samples */
	float dz;		/* z sampling interval */
	float fz;		/* first z sample */
	int live;		/* random number used to denote a live cell */
	int dead;		/* random number used to denote a dead cell */
	float wmin;		/* minimum (reference) frequency */
	float lmin;		/* minimum beamwidth for frequency wmin */
	Cell **cell;	/* cell array[mx][mz] */
	Ray *ray;		/* ray */
} ShotCells;

/* Ray functions */
Ray *makeRay (float x0, float z0, float a0, int nt, float dt, float ft,
	float **a1111xz, float **a3333xz, float **a1133xz, float **a1313xz, 
	float **a1113xz, float **a3313xz,
	int nx, float dx, float fx, int nz, float dz, float fz, int type);
void freeRay (Ray *ray);
int nearestRayStep (Ray *ray, float x, float z);

/* Velocity functions */
void *vel2Alloc (int nx, float dx, float fx,
	int nz, float dz, float fz, float **v);
void vel2Free (void *vel2);
void vel2Interp (void *vel2, float x, float z,
	float *v, float *vx, float *vz, float *vxx, float *vxz, float *vzz);

void accBeam (Ray *ray, float fmin, float lmin,
	int nt, float dt, float ft, float *f,
	int nx, int snx, float dx, float fx, int nz, float dz, float fz, ShotCells *shotcells, float **g);
ShotCells *shotrayTimeAmp (Ray *ray, float fmin, float lmin,
	int nt, float dt, float ft, int nx, float dx, float fx,
	int nz, float dz, float fz);

/* functions defined and used internally */
static void premiggban(float bwh, float fmin, float fmax, float amin, float amax, int nt, float dt,
				int nx, int snx, float sfx, float lx, float rx,  float dx, float dis_trace, int nz, 
				float dz, float **a1111, float **a3333,float **a1133, float **a1313, float **a1113,
				float **a3313, float **f, float **g);

int main()
{
	int i,j,k,m;
	int num;
	int nx;		/* number of traces in input		*/
	int nz; 	/* number of depth values in output	*/
	int nt;		/* number of time samples in input	*/
	int ix,iz;	/* counters				*/
	int lsnx;
	int snx;    /*number of sigle-shot receiver*/
	int n_shot; /*the number of shot*/
	int ia1111=1,ia1133=1,ia1313=1,ia1113=1,ia3313=1,
		idelta=1,iepsilon=1;

	float dx,dz,dt,fmin,fmax,amin,amax,vavg,vavg1,vavg3,bwh,shot_itv,f_shot,sfx,dis_trace,lx,rx, **a1111,**a3333,
		**a1133,**a1313,**a1113,**a3313,**delta,**epsilon,**f,**g;

	int trace_num,samples_per_trace;
	float **trace_data;
	Segy_volumn_hdr volumn_hdr;
	Segy_volumn_hdr write_volumn_hdr;
	Segy_trace_hdr *trace_hdr;
	Segy_trace_hdr *write_trace_hdr;
	int *sinshot_tracenum;      //the trace number of every single shot

	char a1111file[100],a3333file[100],a1133file[100],a1313file[100],
		a1113file[100],a3313file[100],deltafile[100],epsilonfile[100];
	char seifile[100], migfile[100];
	char buff[200];

	FILE *a1111fp,*a3333fp,*a1133fp,
		*a1313fp,*a1113fp,*a3313fp,*deltafp,*epsilonfp;

	FILE *fp=NULL;
	FILE *tracefp;	/* seimic record file */
	FILE *migfp;    /* migration file */

	/*read the canshu file*/
	fp=fopen("canshu.txt","r");
	if(fp==NULL)
	{
		printf("canshu.txt is not exist !\n");
		return ;
	}
    printf("canshu.txt has been opened !\n");
	fgets(buff,200,fp);
	fscanf(fp,"%d",&nx);fscanf(fp,"%d",&nz);fscanf(fp,"%f",&dx);fscanf(fp,"%f",&dz);
	fgets(buff,200,fp);fgets(buff,200,fp);
	fscanf(fp,"%f",&dt);fscanf(fp,"%d",&nt);
	fgets(buff,200,fp);fgets(buff,200,fp);
	fscanf(fp,"%f",&dis_trace);
	fgets(buff,200,fp);fgets(buff,200,fp);
	fscanf(fp,"%d",&n_shot);
	fgets(buff,200,fp);fgets(buff,200,fp);
	fscanf(fp,"%f",&shot_itv);
	fgets(buff,200,fp);fgets(buff,200,fp);
	fscanf(fp,"%f",&f_shot);
    fgets(buff,200,fp);fgets(buff,200,fp);
    fscanf(fp,"%f",&fmin);fscanf(fp,"%f",&fmax);
	fgets(buff,200,fp);fgets(buff,200,fp);
	fscanf(fp,"%s", &a1111file);
	fgets(buff,200,fp);fgets(buff,200,fp);
	fscanf(fp,"%s", &a3333file);
	fgets(buff,200,fp);fgets(buff,200,fp);
	fscanf(fp,"%s", &a1133file);
	fgets(buff,200,fp);fgets(buff,200,fp);
	fscanf(fp,"%s", &a1313file);
	fgets(buff,200,fp);fgets(buff,200,fp);
	fscanf(fp,"%s", &a1113file);
	fgets(buff,200,fp);fgets(buff,200,fp);
	fscanf(fp,"%s", &a3313file);
	fgets(buff,200,fp);fgets(buff,200,fp);
	fscanf(fp,"%s", &epsilonfile);
	fgets(buff,200,fp);fgets(buff,200,fp);
	fscanf(fp,"%s", &deltafile);
	fgets(buff,200,fp);fgets(buff,200,fp);
	fscanf(fp,"%s", &seifile);
	fgets(buff,200,fp);fgets(buff,200,fp);
	fscanf(fp,"%s", &migfile);
	fclose(fp);

	amax = 65.0;
	amin = -amax;

	/* allocate workspace */
	a1111 = alloc2float(nz,nx);
	a3333 = alloc2float(nz,nx);
	a1133 = alloc2float(nz,nx);
	a1313 = alloc2float(nz,nx);
	a1113 = alloc2float(nz,nx);
	a3313 = alloc2float(nz,nx);
	delta = alloc2float(nz,nx);
	epsilon = alloc2float(nz,nx);
	g = alloc2float(nz,nx);

	//--------------- read seismic record(segy data)----------------//
	volumn_hdr=Readvolumn(seifile);

	samples_per_trace=volumn_hdr.samples_per_trace;
	trace_num=get_trace_num(samples_per_trace, seifile);

	/* allocate workspace */
	trace_data = alloc2float(samples_per_trace,trace_num);
	trace_hdr = (Segy_trace_hdr*)alloc1(trace_num,sizeof(Segy_trace_hdr));
	
	/* read segy file */
	Readsegy(trace_num, samples_per_trace, trace_hdr,trace_data, seifile);

	//------------------the trace number of every shot-----------------//
	sinshot_tracenum=alloc1int(n_shot);

	num=0;
	for(i=1;i<trace_num;i++)
	{
		if(trace_hdr[i].trace_sequence_number_within_original_field_record==1)
		{
			sinshot_tracenum[num]=trace_hdr[i-1].trace_sequence_number_within_original_field_record;
			num=num+1;
		}
	}
	sinshot_tracenum[n_shot-1]=trace_hdr[trace_num-1].trace_sequence_number_within_original_field_record;

	/*anisotropic parameters*/
	if ((a1111fp=fopen(a1111file,"rb"))==NULL) 
			ia1111=0;
	if ((a1133fp=fopen(a1133file,"rb"))==NULL) 
			ia1133=0;
	if ((a1313fp=fopen(a1313file,"rb"))==NULL) 
			ia1313=0;
	if ((a1113fp=fopen(a1113file,"rb"))==NULL) 
			ia1113=0;
	if ((a3313fp=fopen(a3313file,"rb"))==NULL) 
			ia3313=0;
	if ((deltafp=fopen(deltafile,"rb"))==NULL) 
			idelta=0;
	if ((epsilonfp=fopen(epsilonfile,"rb"))==NULL) 
			iepsilon=0;

	
	/* read required velocities*/
	a3333fp=fopen(a3333file,"rb");
	for (ix=0; ix<nx; ++ix)
	{
	    fread(a3333[ix],sizeof(float), nz, a3333fp);
	}
	fclose(a3333fp);

	/* read and halve velocities*/
	if(ia1111){
		for (ix=0; ix<nx; ++ix)
			fread(a1111[ix],sizeof(float), nz, a1111fp);
	}else{
		for (ix=0; ix<nx; ++ix)
			for (iz=0; iz<nz; ++iz)
				a1111[ix][iz] = a3333[ix][iz];
	}		
	
	if(ia1133){
		for (ix=0; ix<nx; ++ix)
			fread(a1133[ix],sizeof(float), nz, a1133fp);
	}else{
		for (ix=0; ix<nx; ++ix)
			for (iz=0; iz<nz; ++iz)
				a1133[ix][iz] = .4*a3333[ix][iz];
	}

	if(ia1313){
		for (ix=0; ix<nx; ++ix)
			fread(a1313[ix],sizeof(float), nz, a1313fp);
	}else{
		for (ix=0; ix<nx; ++ix)
			for (iz=0; iz<nz; ++iz)
				a1313[ix][iz] = .3*a3333[ix][iz];
	}

	if(ia1113){
		for (ix=0; ix<nx; ++ix)
			fread(a1113[ix],sizeof(float), nz, a1113fp);
	}else{
		for (ix=0; ix<nx; ++ix)
			for (iz=0; iz<nz; ++iz)
				a1113[ix][iz] = 0;
	}

	if(ia3313){
		for (ix=0; ix<nx; ++ix)
			fread(a3313[ix],sizeof(float), nz, a3313fp);
	}else{
		for (ix=0; ix<nx; ++ix)
			for (iz=0; iz<nz; ++iz)
				a3313[ix][iz] = 0;
	}

	/* if specified read Thomsen's parameters*/
	if(idelta==1 || iepsilon==1) {
		if(idelta){
			for (ix=0; ix<nx; ++ix)
			fread(delta[ix],sizeof(float), nz, deltafp);
		}else{
			for (ix=0; ix<nx; ++ix)
				for (iz=0; iz<nz; ++iz)
					delta[ix][iz] = 0;
		}
		if(iepsilon){
			for (ix=0; ix<nx; ++ix)
			fread(epsilon[ix],sizeof(float), nz, epsilonfp);
		}else{
			for (ix=0; ix<nx; ++ix)
				for (iz=0; iz<nz; ++iz)
					epsilon[ix][iz] = 0;
		}
		for (ix=0; ix<nx; ++ix)
			for (iz=0; iz<nz; ++iz){
				a1111[ix][iz] = a3333[ix][iz]+2*epsilon[ix][iz]
					*a3333[ix][iz];
				a1133[ix][iz] = sqrt(2*delta[ix][iz]*a3333[ix][iz]*
					(a3333[ix][iz]-a1313[ix][iz])+(a3333[ix][iz]-
					a1313[ix][iz])*(a3333[ix][iz]-a1313[ix][iz]))
					-a1313[ix][iz];
			}
	}
	
	/* free workspace */
	free2float(delta);
	free2float(epsilon);
		

	/*determine average velocity*/
	vavg3=0.0; 
	for (ix=0,vavg1=0.0; ix<nx; ++ix) {
		for (iz=0; iz<nz; ++iz) {
			vavg1 += sqrt(a1111[ix][iz]);
			vavg3 += sqrt(a3333[ix][iz]);
		}
	}
	vavg = (vavg1+vavg3)/2;
	vavg /= nx*nz;
	
	/* get beam half-width */
	bwh = vavg/fmin;
	printf("%f\n",bwh);

	/* migrated image */
	for (ix=0; ix<nx; ++ix)
		for (iz=0; iz<nz; ++iz)
			g[ix][iz] = 0.0;

	m=0;
	for(k=0;k<n_shot;k++)
	{
		snx=sinshot_tracenum[k];
		f = alloc2float(nt,snx);
		sfx=f_shot+k*shot_itv;
		lx=ABS(trace_hdr[m].distance_from_source_point_to_receiver_group);
		rx=ABS(trace_hdr[m+snx-1].distance_from_source_point_to_receiver_group);
		lsnx=lx/dis_trace;
		for(i=0;i<snx;i++)
		{
			for(j=0;j<nt;j++)
			{
				f[i][j]=trace_data[m+i][j];
			}
		}

		/* rectify the polarity reversal of converted wave */
		for(i=0;i<lsnx;i++)
		{
			for(j=0;j<nt;j++)
			{
				f[i][j]=-f[i][j];
			}
		}

		fprintf(stderr,"the number of shot record:n_shot=%d\n",k);

		/* migrate */
	    premiggban(bwh,fmin,fmax,amin,amax,nt,dt,nx,snx,sfx,lx,rx,dx,dis_trace,nz,dz,a1111,a3333,
		      a1133,a1313,a1113,a3313,f,g);

		m=m+snx;
		free2float(f);
	}

	//------------write migration data(sgy data)--------//
	write_trace_hdr = (Segy_trace_hdr*)alloc1(nx,sizeof(Segy_trace_hdr));
	/*the data of volumn header*/
	write_volumn_hdr=volumn_hdr;
	write_volumn_hdr.sample_data_interval_ms=dz*1000;
	write_volumn_hdr.original_data_interval_ms=dz*1000;
	write_volumn_hdr.samples_per_trace=nz;
	write_volumn_hdr.original_samples_per_trace=nz;
	/*the data of trace header*/
	write_trace_hdr[0]=trace_hdr[0];
	write_trace_hdr[0].distance_from_source_point_to_receiver_group=0;
	write_trace_hdr[0].x_source_coordinate=0;
	write_trace_hdr[0].samples_in_this_trace=nz;
	write_trace_hdr[0].sample_intervall=dz*1000;
	for(i=1;i<nx;i++)
	{
		write_trace_hdr[i]=write_trace_hdr[0];
		write_trace_hdr[i].trace_sequence_number_within_line=i+1;
	}

	Writesegy(nx, nz, write_volumn_hdr, write_trace_hdr,g, migfile);


	/* free workspace */
	free2float(a1111);
	free2float(a3333);
	free2float(a1133);
	free2float(a1313);
	free2float(a1113);
	free2float(a3313);
	free2float(g);
	free1int(sinshot_tracenum);
	
	return EXIT_SUCCESS;
}

static void premiggban (float bwh, float fmin, float fmax, float amin, float amax, int nt, float dt,
				int nx, int snx, float sfx, float lx, float rx,  float dx, float dis_trace, int nz, 
				float dz, float **a1111, float **a3333,float **a1133, float **a1313, float **a1113,
				float **a3313, float **f, float **g)
/*****************************************************************************
Migrate single shot record data via accumulation of Gaussian beams for anistropic.
******************************************************************************
Input:
bwh		horizontal beam half-width at surface z=0
fmin		minimum frequency (cycles per unit time)
fmax		maximum frequency (cycles per unit time)
amin		minimum emergence angle at surface z=0 (degrees)
amax		maximum emergence angle at surface z=0 (degrees)
nt		number of time samples
dt		time sampling interval (first time assumed to be zero)
nx		number of x samples
dx		x sampling interval
nz		number of z samples
dz		z sampling interval
f		array[nx][nt] containing zero-offset data f(t,x)

*****************************************************************************
Output:
g		array[nx][nz] containing migrated image g(x,z)
*****************************************************************************/
{
	int nxb,npx,ntau,ipx,ix,ixb,ixlo,ixhi,nxw,iz,fxn,mul,ixm,i,bcxn;
	int ir;
	float ft,fx,fz,xwh,dxb,fxb,xb,vmin,dpx,fpx,px,
		taupad,dtau,ftau,fxw,pxmin,pxmax,
		a0,x0,z0,bwhc,**b;

	float f1,f2,f3,f4,f5,f6,alpha,beta,gamma,det,
		rad,signbeta,eps,pz2,q,vp,px2,pz;
	float aa,bb,ccc,dd,ee,fxx,dfx,f7,f8,f9,f10,f11,dpz;
	float gamm11,gamm33,gamm13,vpmin,vpmax,vp2,s,c,ss,cc;

	Ray *ray;
	/* define shot ray variable */
	int jpx,npxh;
	float a01,x01,z01,px1,bwhc1;
	Ray *ray1;
	ShotCells *shotcells;

	/* define the starting point of migration */
	float befx; //the starting point of migration
	int msnx;   //number of x samples in single shot migration area

	/* define left and right of shot extending migration area */
	float l_aper; //extending left migration area
	float r_aper; //extending right migration area
	
	float polx,porx; //the position(coordinate) of first and last receiver
	int polxn;

	float rayfx;    //the origination position of ray trace
	rayfx=0.0;

	mul=dis_trace/dx;

	/* first t, x(m), and z(m)*/
	ft=fz=0.0;
	fx=sfx;  //the shot posizion coordinate
	fxn=fx/dx;

	/* the position of first and last receiver */
	polx=fx-lx;
	porx=fx+rx;
	polxn=polx/dx;
	
	/* present extending migration area and calculate the starting point of migration,
	   extending in left and right of shot area, the size extending depending on model. */
	l_aper=300;
	if(polx-l_aper<0) l_aper=polx;
	r_aper=300;
	if(porx+r_aper>(nx-1)*dx) r_aper=(nx-1)*dx-porx;
	befx=polx-l_aper;
	if(befx<0) befx=0.0;
	msnx=snx*dis_trace/dx-1+l_aper/dx+r_aper/dx;
	if(msnx>nx) msnx=nx;

	/* convert minimum and maximum angles to radians */
	amin *= PI/180.0;
	amax *= PI/180.0;
	if (amin>amax) {
		float atemp=amin;
		amin = amax;
		amax = atemp;
	}
	
	/* minimum velocity at surface z=0 */
	ixm = 0;
	for (ix=1,vmin=a3333[0][0]; ix<nx; ++ix)
		if (a3333[ix][0]<vmin) { 
			ixm=ix;
			vmin = a3333[ix][0];
		}

	/* phase velocity at min. and max. angles */
	c = cos(amin);
	s = sin(amin);
	cc = c*c;
	ss = s*s;
	gamm11 = a1111[ixm][0]*ss+ a1313[ixm][0]*cc +2*a1113[ixm][0]*s*c;
	gamm33 = a3333[ixm][0]*cc + a1313[ixm][0]*ss+2*a3313[ixm][0]*s*c;
	gamm13 = (a1133[ixm][0]+a1313[ixm][0])*s*c+ a1113[ixm][0]*ss+ 
		a3313[ixm][0]*cc;
	vp2    = gamm11+gamm33+sqrt((gamm11+gamm33)*(gamm11+gamm33)-
			4*(gamm11*gamm33-gamm13*gamm13));
	vpmin     = sqrt(vp2*.5);

	c = cos(amax);
	s = sin(amax);
	cc = c*c;
	ss = s*s;
	gamm11 = a1111[ixm][0]*ss+ a1313[ixm][0]*cc +2*a1113[ixm][0]*s*c;
	gamm33 = a3333[ixm][0]*cc + a1313[ixm][0]*ss+2*a3313[ixm][0]*s*c;
	gamm13 = (a1133[ixm][0]+a1313[ixm][0])*s*c+ a1113[ixm][0]*ss+ 
		a3313[ixm][0]*cc;
	vp2    = gamm11+gamm33+sqrt((gamm11+gamm33)*(gamm11+gamm33)-
			4*(gamm11*gamm33-gamm13*gamm13));
	vpmax     = sqrt(vp2*.5);

	/* beam sampling */
	pxmin = sin(amin)/vpmin;
	pxmax = sin(amax)/vpmax;
	dpx = sqrt(vmin)*1.0/(6.0*bwh*sqrt(fmin*fmax)*.5*(vpmin+vpmax));
	npx = 1+(pxmax-pxmin)/dpx;
	fpx = pxmin+0.5*(pxmax-pxmin-(npx-1)*dpx);


	fprintf(stderr,"npx=%d dpx=%g fpx=%g\n",npx,dpx,fpx);

	/* loop over reviewer*/
	for (ir=0; ir<snx; ++ir) {

		fprintf(stderr,"ir = %d\n",ir);

		/* loop over shot point */
		for (ipx=0,px1=fpx;ipx<npx;++ipx,px1+=dpx){

			/*calculating emergence angle from ray parameter*/
			f1 = a1111[fxn][0]+a1313[fxn][0];
			f2 = a3333[fxn][0]+a1313[fxn][0];
			f3 = a1111[fxn][0]-a1313[fxn][0];
			f4 = a3333[fxn][0]-a1313[fxn][0];
			f5 = 2.*(a1133[fxn][0]+a1313[fxn][0])*
				(a1133[fxn][0]+a1313[fxn][0]);
			f6 = 2;

			eps   = .0001;
			px2   = px1*px1;
	    		alpha = f2*f2-f4*f4;
	    		beta  = 2*((f1*f2+f3*f4-f5)*px2-f2*f6);
	    		gamma = f6*f6-(2.*f1*f6-(f1*f1-f3*f3)*px2)*px2;
	    		det   = beta*beta-4.*alpha*gamma;

			if (det<0) continue;
			rad = sqrt(det);
			if(ABS(beta)>eps)   signbeta = ABS(beta)/beta;
	    		else                signbeta = 1.;
			q    = -.5*(beta+signbeta*rad);
			pz2  = gamma/q;
			if(pz2<0) continue;
			pz   = sqrt(pz2);

			f7   = 2*(a1113[fxn][0]+a3313[fxn][0]);
			f8   = 2*(a1113[fxn][0]+a3313[fxn][0]);		
			f9   = a3313[fxn][0];
			f10  = a1113[fxn][0];
			f11  = 2*(a1133[fxn][0]+a1313[fxn][0]);
			aa   = f2*f2-f4*f4-f9*f9;
			bb   = (2*f2*f7+2*f4*f8-2*f11*f9)*px1;
			ccc  = f7*f7*px2-2*(f6-f1*px2)*f2-f8*f8*px2+
				2*f3*f4*px2-f11*f11*px2-2*f10*f9*px2;
			dd   = -2*(f6-f1*px2)*f7*px1-2*f3*f8*px2*px1-
				2*f11*f10*px1*px2;
			ee   = (f6-f1*px2)*(f6-f1*px2)-f3*f3*px2*px2-
				f10*f10*px2*px2;
			for (i=0; i<20; i++) {
				fxx = aa*pz2*pz2+bb*pz2*pz+ccc*pz2+
					dd*pz+ee;
				dfx = 4*aa*pz*pz2+3*bb*pz2+2*ccc*pz+dd;
				if(dfx==0) break;
				dpz = -fxx/dfx;
				if(ABS(dpz)<.000001) break;
				pz += dpz;
				pz2 = pz*pz;
			}
			if (pz<0) continue;
			vp   = 1/sqrt(px2+pz2);

			/* sine of emergence angle; skip if out of bounds */
			if (px1*vp>sin(amax)+0.01) continue;
			if (px1*vp<sin(amin)-0.01) continue;

			/* emergence angle and location */
			a01 = asin(px1*vp);
			//a01 = asin(px1*vpmin);
			x01 = fx;
			z01 = fz;

			/* beam half-width adjusted for cosine of angle */
			bwhc1 = bwh*cos(a01);

			/* trace ray */
			ray1 = makeRay(x01,z01,a01,nt,dt,ft,
				a1111,a3333,a1133,a1313,a1113,a3313,
				nx,dx,rayfx,nz,dz,fz,0);
			
			/* accumulate its Gaussian beam complex time and amplitude for source ray */
			shotcells = shotrayTimeAmp (ray1,fmin,bwhc1,nt,dt,ft,msnx,dx,befx,nz,dz,fz);

			/* loop over beams */
			for (jpx=0,px=fpx; jpx<npx; ++jpx,px+=dpx) {

				bcxn=ir*mul+polxn;   //beam centre position
				/*calculating emergence angle from ray parameter*/
				f1 = a1111[bcxn][0]+a1313[bcxn][0];
				f2 = a3333[bcxn][0]+a1313[bcxn][0];
				f3 = a1111[bcxn][0]-a1313[bcxn][0];
				f4 = a3333[bcxn][0]-a1313[bcxn][0];
				f5 = 2.*(a1133[bcxn][0]+a1313[bcxn][0])*
					(a1133[bcxn][0]+a1313[bcxn][0]);
				f6 = 2;

				eps   = .0001;
				px2   = px*px;
	    			alpha = f2*f2-f4*f4;
	    			beta  = 2*((f1*f2+f3*f4-f5)*px2-f2*f6);
	    			gamma = f6*f6-(2.*f1*f6-(f1*f1-f3*f3)*px2)*px2;
	    			det   = beta*beta-4.*alpha*gamma;

				if (det<0) continue;
				rad = sqrt(det);
				if(ABS(beta)>eps)   signbeta = ABS(beta)/beta;
	    			else                signbeta = 1.;
				q    = -.5*(beta-signbeta*rad);
				pz2  = gamma/q;
				if(pz2<0) continue;
				pz   = sqrt(pz2);

				f7   = 2*(a1113[bcxn][0]+a3313[bcxn][0]);
				f8   = 2*(a1113[bcxn][0]+a3313[bcxn][0]);		
				f9   = a3313[bcxn][0];
				f10  = a1113[bcxn][0];
				f11  = 2*(a1133[bcxn][0]+a1313[bcxn][0]);
				aa   = f2*f2-f4*f4-f9*f9;
				bb   = (2*f2*f7+2*f4*f8-2*f11*f9)*px;
				ccc  = f7*f7*px2-2*(f6-f1*px2)*f2-f8*f8*px2+
					2*f3*f4*px2-f11*f11*px2-2*f10*f9*px2;
				dd   = -2*(f6-f1*px2)*f7*px-2*f3*f8*px2*px-
					2*f11*f10*px*px2;
				ee   = (f6-f1*px2)*(f6-f1*px2)-f3*f3*px2*px2-
					f10*f10*px2*px2;
				for (i=0; i<20; i++) {
					fxx = aa*pz2*pz2+bb*pz2*pz+ccc*pz2+
						dd*pz+ee;
					dfx = 4*aa*pz*pz2+3*bb*pz2+2*ccc*pz+dd;
					if(dfx==0) break;
					dpz = -fxx/dfx;
					if(ABS(dpz)<.000001) break;
					pz += dpz;
					pz2 = pz*pz;
				}
				if (pz<0) continue;
				vp   = 1/sqrt(px2+pz2);

				/* sine of emergence angle; skip if out of bounds */
				if (px*vp>sin(amax)+0.01) continue;
				if (px*vp<sin(amin)-0.01) continue;

				/* emergence angle and location */
				a0 = -asin(px*vp);
				//a0 = -asin(px*vpmin);
				x0 = polx+ir*dis_trace;
				z0 = fz;

				/* beam half-width adjusted for cosine of angle */
				bwhc = bwh*cos(a0);

				/* trace ray */
				ray = makeRay(x0,z0,a0,nt,dt,ft,
				a1111,a3333,a1133,a1313,a1113,a3313,
				nx,dx,rayfx,nz,dz,fz,1);

				/* accumulate contribution of beam in migrated image */
				accBeam(ray,fmin,bwhc,
					nt,dt,ft,f[ir],
					nx,msnx,dx,befx,nz,dz,fz,shotcells,g);
				
				/* free ray */
				freeRay(ray);
			}

			freeRay(ray1);
			/* free shotcells */
			free2((void**)shotcells->cell);
			free1((void*)shotcells);
		}
	}
}

/* circle for efficiently finding nearest ray step */
typedef struct CircleStruct {
	int irsf;               /* index of first raystep in circle */
	int irsl;               /* index of last raystep in circle */
	float x;                /* x coordinate of center of circle */
	float z;                /* z coordinate of center of circle */
	float r;                /* radius of circle */
} Circle;

/* functions defined and used internally */
static Circle *makeCircles (int nc, int nrs, RayStep *rs);

int nearestRayStep (Ray *ray, float x, float z)
/*****************************************************************************
Determine index of ray step nearest to point (x,z).
******************************************************************************
Input:
ray		ray
x		x coordinate
z		z coordinate
*****************************************************************************
Returned:	index of nearest ray step
*****************************************************************************
Credits: CWP: Dave Hale
*****************************************************************************/
{
	int nrs=ray->nrs,ic=ray->ic,nc=ray->nc;
	RayStep *rs=ray->rs;
	Circle *c=ray->c;
	int irs,irsf,irsl,irsmin=0,update,jc,js,kc;
	float dsmin,ds,dx,dz,dmin,rdmin,xrs,zrs;

	/* if necessary, make circles localizing ray steps */
	if (c==NULL) {
		ray->ic = ic = 0;
		ray->nc = nc = sqrt((float)nrs);
		ray->c = c = makeCircles(nc,nrs,rs);
	}
	
	/* initialize minimum distance and minimum distance-squared */
	dx = x-c[ic].x;
	dz = z-c[ic].z;
	dmin = 2.0*(sqrt(dx*dx+dz*dz)+c[ic].r);
	dsmin = dmin*dmin;

	/* loop over all circles */
	for (kc=0,jc=ic,js=0; kc<nc; ++kc) {
		
		/* distance-squared to center of circle */
		dx = x-c[jc].x;
		dz = z-c[jc].z;
		ds = dx*dx+dz*dz;

		/* radius of circle plus minimum distance (so far) */
		rdmin = c[jc].r+dmin;

		/* if circle could possible contain a nearer ray step */
		if (ds<=rdmin*rdmin) {

			/* search circle for nearest ray step */ 
			irsf = c[jc].irsf;
			irsl = c[jc].irsl;
			update = 0;
			for (irs=irsf; irs<=irsl; ++irs) {
				xrs = rs[irs].x;
				zrs = rs[irs].z;
				dx = x-xrs;
				dz = z-zrs;
				ds = dx*dx+dz*dz;
				if (ds<dsmin) {
					dsmin = ds;
					irsmin = irs;
					update = 1;
				}
			}

			/* if a nearer ray step was found inside circle */
			if (update) {

				/* update minimum distance */
				dmin = sqrt(dsmin);

				/* remember the circle */
				ic = jc;
			}
		}
		
		/* search circles in alternating directions */
		js = (js>0)?-js-1:-js+1;
		jc += js;
		if (jc<0 || jc>=nc) {
			js = (js>0)?-js-1:-js+1;
			jc += js;
		}
	}

	/* remember the circle containing the nearest ray step */
	ray->ic = ic;

	if (irsmin<0 || irsmin>=nrs)
		fprintf(stderr,"irsmin=%d\n",irsmin);

	/* return index of nearest ray step */
	return irsmin;
}


Circle *makeCircles (int nc, int nrs, RayStep *rs)
/*****************************************************************************
Make circles used to speed up determination of nearest ray step.
******************************************************************************
Input:
nc		number of circles to make
nrs		number of ray steps
rs		array[nrs] of ray steps
*****************************************************************************
Returned:	array[nc] of circles
*****************************************************************************/
{
	int nrsc,ic,irsf,irsl,irs;
	float xmin,xmax,zmin,zmax,x,z,r;
	Circle *c;

	/* allocate space for circles */
	c = (Circle*)alloc1(nc,sizeof(Circle));

	/* determine typical number of ray steps per circle */
	nrsc = 1+(nrs-1)/nc;

	/* loop over circles */
	for (ic=0; ic<nc; ++ic) {
		
		/* index of first and last raystep */
		irsf = ic*nrsc;
		irsl = irsf+nrsc-1;
		if (irsf>=nrs) irsf = nrs-1;
		if (irsl>=nrs) irsl = nrs-1;

		/* coordinate bounds of ray steps */
		xmin = xmax = rs[irsf].x;
		zmin = zmax = rs[irsf].z;
		for (irs=irsf+1; irs<=irsl; ++irs) {
			if (rs[irs].x<xmin) xmin = rs[irs].x;
			if (rs[irs].x>xmax) xmax = rs[irs].x;
			if (rs[irs].z<zmin) zmin = rs[irs].z;
			if (rs[irs].z>zmax) zmax = rs[irs].z;
		}

		/* center and radius of circle */
		x = 0.5*(xmin+xmax);
		z = 0.5*(zmin+zmax);
		r = sqrt((x-xmin)*(x-xmin)+(z-zmin)*(z-zmin));

		/* set circle */
		c[ic].irsf = irsf;
		c[ic].irsl = irsl;
		c[ic].x = x;
		c[ic].z = z;
		c[ic].r = r;
	}

	return c;
}

/************************************************************************
* Functions for Gaussian Beam computation		                *
*************************************************************************/

/* size of cells in which to linearly interpolate complex time and amplitude */
#define CELLSIZE 8

/* factor by which to oversample time for linear interpolation of traces */
#define NOVERSAMPLE 4

/* number of exponential decay filters */
#define NFILTER 6

/* exp(EXPMIN) is assumed to be negligible */
#define EXPMIN (-5.0)

/* functions defined and used internally */
static BeamData* beamData (float wmin, int nt, float dt, float ft, float *f);
static void shotsetCell(ShotCells *shotcells, int jx,int jz);
static void setCell (Cells *cells, int jx, int jz);
static int shotcellTimeAmp (ShotCells *shotcells, int jx,int jz);
static int cellTimeAmp (Cells *cells, int jx, int jz);

static void accCell (Cells *cells, ShotCells *shotcells);
static void cellBeam (Cells *cells, ShotCells *shotcells, int jx, int jz);


ShotCells *shotrayTimeAmp (Ray *ray, float fmin, float lmin,
	int nt, float dt, float ft, int nx, float dx, float fx,
	int nz, float dz, float fz)
/*********************************************************************************
Set a cell by computing its Gaussian beam complex time and amplitude in source ray.  
**********************************************************************************
Input:
ray		ray parameters sampled at discrete ray steps
fmin	minimum frequency (cycles per unit time)
lmin	initial beam width for frequency wmin
nt		number of time samples
dt		time sampling interval
ft		first time sample
nx		number of x samples in single shot migration area
dx		x sampling interval
fx		first x sample
nz		number of z samples
dz		z sampling interval
fz		first z sample
*****************************************************************************/
{
	int lx,lz,mx,mz,jx,jz,live,dead;
	float wmin;
	RayStep *rs=ray->rs;
	ShotCells *shotcells;
	Cell **cell;

	/* frequency in radians per unit time */
	wmin = 2.0*PI*fmin;

	/* random numbers used to denote live and dead cells */
	live = 1+(int)(1.0e7*franuni());
	dead = 1+(int)(1.0e7*franuni());
	
	/* number of samples per cell */
	lx = CELLSIZE;
	lz = CELLSIZE;

	/* number of cells */
	mx = 2+(nx-1)/lx;
	mz = 2+(nz-1)/lz;

	/* allocate cells */
	shotcells = (ShotCells*)alloc1(1,sizeof(ShotCells));
	cell = (Cell**)alloc2(mz,mx,sizeof(Cell));

	/* set information needed to set and fill cells */
	shotcells->nt = nt;
	shotcells->dt = dt;
	shotcells->ft = ft;
	shotcells->lx = lx;
	shotcells->mx = mx;
	shotcells->nx = nx;
	shotcells->dx = dx;
	shotcells->fx = fx;
	shotcells->lz = lz;
	shotcells->mz = mz;
	shotcells->nz = nz;
	shotcells->dz = dz;
	shotcells->fz = fz;
	shotcells->live = live;
	shotcells->dead = dead;
	shotcells->wmin = wmin;
	shotcells->lmin = lmin;
	shotcells->cell = cell;
	shotcells->ray = ray;
	
	/* cell closest to initial point on ray will be first live cell */
	jx = NINT((rs[0].x-fx)/dx/lx);
	jz = NINT((rs[0].z-fz)/dz/lz);
	
	/* set first live cell and its neighbors recursively */
	shotsetCell(shotcells,jx,jz);

	return shotcells;
}

static void shotsetCell (ShotCells *shotcells, int jx, int jz)
/************************************************************************************
Set a cell by computing its Gaussian beam complex time and amplitude in shot position.  
If the amplitude is non-zero, set neighboring cells recursively.
*************************************************************************************
Input:
shotcells		pointer to cells
jx				x index of the cell to set
jz				z index of the cell to set
******************************************************************************
Notes:
To reduce the amount of memory required for recursion, the actual 
computation of complex time and amplitude is performed by the shotcellTimeAmp() 
function, so that no local variables are required in this function, except
for the input arguments themselves.
*****************************************************************************/
{
	/* if cell is out of bounds, return */
	if (jx<0 || jx>=shotcells->mx || jz<0 || jz>=shotcells->mz) return;

	/* if cell is live, return */
	if (shotcells->cell[jx][jz].live==shotcells->live) return;

	/* make cell live */
	shotcells->cell[jx][jz].live = shotcells->live;

	/* compute complex time and amplitude.  If amplitude is
	 * big enough, recursively set neighboring cells. */
	if (shotcellTimeAmp(shotcells,jx,jz)) {
		shotsetCell(shotcells,jx+1,jz);
		shotsetCell(shotcells,jx-1,jz);
		shotsetCell(shotcells,jx,jz+1);
		shotsetCell(shotcells,jx,jz-1);
	}
}

static int shotcellTimeAmp (ShotCells *shotcells, int jx, int jz)
/*****************************************************************************
Compute complex time and amplitude for a cell for shot ray.
******************************************************************************
Input:
shotcells		pointer to cells
jx				x index of the cell to set
jz				z index of the cell to set

Returned:	1 if Gaussian amplitude is significant, 0 otherwise
*****************************************************************************/
{
	int lx=shotcells->lx,lz=shotcells->lz;
	float dx=shotcells->dx,fx=shotcells->fx,dz=shotcells->dz,fz=shotcells->fz,
		wmin=shotcells->wmin,lmin=shotcells->lmin;
	Ray *ray=shotcells->ray;
	Cell **cell=shotcells->cell;
	int irs,kmah;
	float tmax,xc,zc,t0,x0,z0,x,z,tr,ti,ar,ai,
		v,q1,p1,q2,p2,e,es,scale,phase,vvmr,vvmi,c,s,ov,
		px,pz,pzpz,pxpx,pxpz,dvdx,dvdz,dvds,dvdn,
		wxxr,wxzr,wzzr,wxxi,wxzi,wzzi;
	RayStep *rs;

	/* maximum time */
	tmax = ray->rs[ray->nrs-1].t;

	/* cell coordinates */
	xc = fx+jx*lx*dx;
	zc = fz+jz*lz*dz;

	/* ray step nearest to cell */
	irs = nearestRayStep(ray,xc,zc);
	rs = &(ray->rs[irs]);
		
	/* ray time and coordinates */
	t0 = rs->t;
	x0 = rs->x;
	z0 = rs->z;
	
	/* real and imaginary parts of v*v*m = v*v*p/q */
	v = rs->v;
	q1 = rs->q1;
	p1 = rs->p1;
	q2 = rs->q2;
	p2 = rs->p2;
	e = wmin*lmin*lmin;
	es = e*e;
	scale = v*v/(q2*q2+es*q1*q1);
	vvmr = (p2*q2+es*p1*q1)*scale;
	vvmi = -e*scale;
	
	/* components of slowness vector, px and pz */
	c = rs->c;
	s = rs->s;
	ov = 1.0/v;
	px = s*ov;
	pz = c*ov;
	pzpz = pz*pz;
	pxpx = px*px;
	pxpz = px*pz;
	
	/* velocity derivatives along tangent and normal */
	dvdx = rs->dvdx;
	dvdz = rs->dvdz;
	dvds = s*dvdx+c*dvdz;
	dvdn = c*dvdx-s*dvdz;
	
	/* real part of W matrix */
	wxxr = pzpz*vvmr-2.0*pxpz*dvdn-pxpx*dvds;
	wxzr = (pxpx-pzpz)*dvdn-pxpz*(vvmr+dvds);
	wzzr = pxpx*vvmr+2.0*pxpz*dvdn-pzpz*dvds;
	
	/* imaginary part of W matrix */
	wxxi = pzpz*vvmi;
	wxzi = -pxpz*vvmi;
	wzzi = pxpx*vvmi;
	
	/* vector from ray to cell */
	x = xc-x0;
	z = zc-z0;
	
	/* real and imaginary parts of complex time */
	tr = t0+(px+0.5*wxxr*x)*x+(pz+wxzr*x+0.5*wzzr*z)*z;
	ti = (0.5*wxxi*x*x+(wxzi*x+0.5*wzzi*z)*z)*2;
	
	/* real and imaginary parts of complex amplitude */
	kmah = rs->kmah;
	scale = pow(es/(q2*q2+es*q1*q1),0.25);
	phase = 0.5*(atan2(q2,q1*e)+2.0*PI*((kmah+1)/2));
	ar = scale*cos(phase);
	ai = scale*sin(phase);

	/* set cell parameters */
	cell[jx][jz].tr = tr;
	cell[jx][jz].ti = ti;
	cell[jx][jz].ar = ar;
	cell[jx][jz].ai = ai;

	/* return 1 if Gaussian amplitude is significant, 0 otherwise */
	return (wmin*ti>EXPMIN && tr<=tmax)?1:0;
}

void accBeam (Ray *ray, float fmin, float lmin,
	int nt, float dt, float ft, float *f,
	int nx, int snx, float dx, float fx, int nz, float dz, float fz, ShotCells *shotcells, float **g)
/*****************************************************************************
Accumulate contribution of one Gaussian beam and source ray.
******************************************************************************
Input:
ray			ray parameters sampled at discrete ray steps
fmin		minimum frequency (cycles per unit time)
lmin		initial beam width for frequency wmin
nt			number of time samples
dt			time sampling interval
ft			first time sample
f			array[nt] containing data for one ray f(t)		
nx			number of x samples in velocity model
snx			number of x samples in single shot migration area
dx			x sampling interval
fx			first x sample
nz			number of z samples
dz			z sampling interval
fz			first z sample
shotcells	pointer to shotcells(complex time and amplitude in shot position)
g			array[nx][nz] in which to accumulate beam

Output:
g		array[nx][nz] after accumulating beam
*****************************************************************************/
{
	int lx,lz,mx,mz,jx,jz,live,dead;
	float wmin;
	RayStep *rs=ray->rs;
	Cells *cells;
	Cell **cell;
	BeamData *bd;

	/* frequency in radians per unit time */
	wmin = 2.0*PI*fmin;

	/* random numbers used to denote live and dead cells */
	live = 1+(int)(1.0e7*franuni());
	dead = 1+(int)(1.0e7*franuni());
	
	/* number of samples per cell */
	lx = CELLSIZE;
	lz = CELLSIZE;

	/* number of cells */
	mx = 2+(snx-1)/lx;
	mz = 2+(nz-1)/lz;

	/* compute complex beam data */
	bd = beamData(wmin,nt,dt,ft,f);

	/* allocate cells */
	cells = (Cells*)alloc1(1,sizeof(Cells));
	cell = (Cell**)alloc2(mz,mx,sizeof(Cell));

	/* set information needed to set and fill cells */
	cells->nt = nt;
	cells->dt = dt;
	cells->ft = ft;
	cells->lx = lx;
	cells->mx = mx;
	cells->nx = nx;
	cells->snx = snx;
	cells->dx = dx;
	cells->fx = fx;
	cells->lz = lz;
	cells->mz = mz;
	cells->nz = nz;
	cells->dz = dz;
	cells->fz = fz;
	cells->live = live;
	cells->dead = dead;
	cells->wmin = wmin;
	cells->lmin = lmin;
	cells->cell = cell;
	cells->ray = ray;
	cells->bd = bd;
	cells->g = g;

	/* cell closest to initial point on ray will be first live cell */
	jx = NINT((rs[0].x-fx)/dx/lx);
	jz = NINT((rs[0].z-fz)/dz/lz);
	
	/* set first live cell and its neighbors recursively */
	setCell(cells,jx,jz);

	/* accumulate beam in first live cell and its neighbors recursively */
	accCell(cells,shotcells);

	/* free complex beam data */
	free2complex(bd->cf);
	free1((void*)bd);

	/* free cells */
	free2((void**)cells->cell);
	free1((void*)cells);
}


static BeamData* beamData (float wmin, int nt, float dt, float ft, float *f)
/*****************************************************************************
Compute filtered complex beam data as a function of real and imaginary time.
******************************************************************************
Input:
wmin		minimum frequency (in radians per unit time)
nt		number of time samples
dt		time sampling interval
ft		first time sample
f		array[nt] containing data to be filtered
*****************************************************************************
Returned:	pointer to beam data
*****************************************************************************/
{
	int ntpad,ntfft,nw,iwnyq,ntrfft,ntr,nti,nwr,it,itr,iti,iw;
	float dw,fw,dtr,ftr,dti,fti,w,ti,scale,*fa;
	complex *ca,*cb,*cfi,**cf;
	BeamData *bd;

	/* pad to avoid wraparound in Hilbert transform */
	ntpad = 25;

	/* fft sampling */
	ntfft = npfaro(nt+ntpad,2*(nt+ntpad));
	nw = ntfft/2+1;
	dw = 2.0*PI/(ntfft*dt);
	fw = 0.0;
	iwnyq = nw-1;

	/* real time sampling (oversample for future linear interpolation) */
	ntrfft = nwr = npfao(NOVERSAMPLE*ntfft,NOVERSAMPLE*ntfft+ntfft);
	dtr = dt*ntfft/ntrfft;
	ftr = ft;
	ntr = 1+(nt+ntpad-1)*dt/dtr;

	/* imaginary time sampling (exponential decay filters) */
	nti = NFILTER;
	dti = EXPMIN/(wmin*(nti-1));
	fti = 0.0;

	/* allocate space for filtered data */
	cf = alloc2complex(ntr,nti);

	/* allocate workspace */
	fa = alloc1float(ntfft);
	ca = alloc1complex(nw);
	cb = alloc1complex(ntrfft);

	/* pad data with zeros */
	for (it=0; it<nt; ++it)
		fa[it] = f[it];
	for (it=nt; it<ntfft; ++it)
		fa[it] = 0.0;

	/* Fourier transform and scale to make complex analytic signal */
	pfarc(1,ntfft,fa,ca);
	for (iw=1; iw<iwnyq; ++iw) {
		ca[iw].r *= 2.0;
		ca[iw].i *= 2.0;
	}

	/* loop over imaginary time */
	for (iti=0,ti=fti; iti<nti; ++iti,ti+=dti) {

		/* apply exponential decay filter */
		for (iw=0,w=fw; iw<nw; ++iw,w+=dw) {
			scale = exp(w*ti);
			cb[iw].r = ca[iw].r*scale;
			cb[iw].i = ca[iw].i*scale;
		}

		/* pad with zeros */
		for (iw=nw; iw<nwr; ++iw)
			cb[iw].r = cb[iw].i = 0.0;

		/* inverse Fourier transform and scale */
		pfacc(-1,ntrfft,cb);
		cfi = cf[iti];
		scale = 1.0/ntfft;
		for (itr=0; itr<ntr; ++itr) {
			cfi[itr].r = scale*cb[itr].r;
			cfi[itr].i = scale*cb[itr].i;
		}
	}

	/* free workspace */
	free1float(fa);
	free1complex(ca);
	free1complex(cb);

	/* return beam data */
	bd = (BeamData*)alloc1(1,sizeof(BeamData));
	bd->ntr = ntr;
	bd->dtr = dtr;
	bd->ftr = ftr;
	bd->nti = nti;
	bd->dti = dti;
	bd->fti = fti;
	bd->cf = cf;
	return bd;
}

static void setCell (Cells *cells, int jx, int jz)
/***********************************************************************************
Set a cell by computing its Gaussian beam complex time and amplitude in beams center.  
If the amplitude is non-zero, set neighboring cells recursively.
************************************************************************************
Input:
cells	pointer to cells
jx		x index of the cell to set
jz		z index of the cell to set
************************************************************************************
Notes:
To reduce the amount of memory required for recursion, the actual 
computation of complex time and amplitude is performed by the cellTimeAmp() 
function, so that no local variables are required in this function, except
for the input arguments themselves.
************************************************************************************/
{
	/* if cell is out of bounds, return */
	if (jx<0 || jx>=cells->mx || jz<0 || jz>=cells->mz) return;

	/* if cell is live, return */
	if (cells->cell[jx][jz].live==cells->live) return;

	/* make cell live */
	cells->cell[jx][jz].live = cells->live;

	/* compute complex time and amplitude.  If amplitude is
	 * big enough, recursively set neighboring cells. */
	if (cellTimeAmp(cells,jx,jz)) {
		setCell(cells,jx+1,jz);
		setCell(cells,jx-1,jz);
		setCell(cells,jx,jz+1);
		setCell(cells,jx,jz-1);
	}
}

static int cellTimeAmp (Cells *cells, int jx, int jz)
/*****************************************************************************
Compute complex and time and amplitude for a cell in beams center.
******************************************************************************
Input:
cells	pointer to cells
jx		x index of the cell to set
jz		z index of the cell to set

Returned:	1 if Gaussian amplitude is significant, 0 otherwise
*****************************************************************************/
{
	int lx=cells->lx,lz=cells->lz;
	float dx=cells->dx,fx=cells->fx,dz=cells->dz,fz=cells->fz,
		wmin=cells->wmin,lmin=cells->lmin;
	Ray *ray=cells->ray;
	Cell **cell=cells->cell;
	int irs,kmah;
	float tmax,xc,zc,t0,x0,z0,x,z,tr,ti,ar,ai,
		v,q1,p1,q2,p2,e,es,scale,phase,vvmr,vvmi,c,s,ov,
		px,pz,pzpz,pxpx,pxpz,dvdx,dvdz,dvds,dvdn,
		wxxr,wxzr,wzzr,wxxi,wxzi,wzzi;
	RayStep *rs;

	/* maximum time */
	tmax = ray->rs[ray->nrs-1].t;

	/* cell coordinates */
	xc = fx+jx*lx*dx;
	zc = fz+jz*lz*dz;

	/* ray step nearest to cell */
	irs = nearestRayStep(ray,xc,zc);
	rs = &(ray->rs[irs]);
		
	/* ray time and coordinates */
	t0 = rs->t;
	x0 = rs->x;
	z0 = rs->z;
	
	/* real and imaginary parts of v*v*m = v*v*p/q */
	v = rs->v;
	q1 = rs->q1;
	p1 = rs->p1;
	q2 = rs->q2;
	p2 = rs->p2;
	e = wmin*lmin*lmin;
	es = e*e;
	scale = v*v/(q2*q2+es*q1*q1);
	vvmr = (p2*q2+es*p1*q1)*scale;
	vvmi = -e*scale;
	
	/* components of slowness vector, px and pz */
	c = rs->c;
	s = rs->s;
	ov = 1.0/v;
	px = s*ov;
	pz = c*ov;
	pzpz = pz*pz;
	pxpx = px*px;
	pxpz = px*pz;
	
	/* velocity derivatives along tangent and normal */
	dvdx = rs->dvdx;
	dvdz = rs->dvdz;
	dvds = s*dvdx+c*dvdz;
	dvdn = c*dvdx-s*dvdz;
	
	/* real part of W matrix */
	wxxr = pzpz*vvmr-2.0*pxpz*dvdn-pxpx*dvds;
	wxzr = (pxpx-pzpz)*dvdn-pxpz*(vvmr+dvds);
	wzzr = pxpx*vvmr+2.0*pxpz*dvdn-pzpz*dvds;
	
	/* imaginary part of W matrix */
	wxxi = pzpz*vvmi;
	wxzi = -pxpz*vvmi;
	wzzi = pxpx*vvmi;
	
	/* vector from ray to cell */
	x = xc-x0;
	z = zc-z0;
	
	/* real and imaginary parts of complex time */
	tr = t0+(px+0.5*wxxr*x)*x+(pz+wxzr*x+0.5*wzzr*z)*z;
	ti = (0.5*wxxi*x*x+(wxzi*x+0.5*wzzi*z)*z)*2;
	
	/* real and imaginary parts of complex amplitude */
	kmah = rs->kmah;
	scale = pow(es/(q2*q2+es*q1*q1),0.25);
	phase = 0.5*(atan2(q2,q1*e)+2.0*PI*((kmah+1)/2));
	ar = scale*cos(phase);
	ai = scale*sin(phase);

	/* set cell parameters */
	cell[jx][jz].tr = tr;
	cell[jx][jz].ti = ti;
	cell[jx][jz].ar = ar;
	cell[jx][jz].ai = ai;

	/* return 1 if Gaussian amplitude is significant, 0 otherwise */
	return (wmin*ti>EXPMIN && tr<=tmax)?1:0;
}

static void accCell (Cells *cells, ShotCells *shotcells)
/*****************************************************************************
Accumulate the contribution of a Gaussian beam in a cell and its
neighbors, recursively.
******************************************************************************
Input:
cells		pointer to cells
jx		x index of the cell to fill
jz		z index of the cell to fill
******************************************************************************
Notes:
To reduce the amount of memory required for recursion, the actual
accumulation is performed by function cellBeam(), so that no local
variables are required in this function, except for the input
arguments themselves.
*****************************************************************************/
{
	int ix,iz,mx,mz;
	mx=cells->mx;
	mz=cells->mz;
	
	for(ix=0;ix<mx;ix++)
	{
		for(iz=0;iz<mz;++iz)
		{
			/* if cell is out of bounds, return */
			if (ix<0 || ix>=cells->mx-1 || iz<0 || iz>=cells->mz-1) break;

			/* if remaining three corners of cell are live */
			if (cells->cell[ix][iz].live==cells->live &&
				cells->cell[ix+1][iz].live==cells->live &&
				cells->cell[ix][iz+1].live==cells->live &&
				cells->cell[ix+1][iz+1].live==cells->live &&
				shotcells->cell[ix][iz].live==shotcells->live &&
				shotcells->cell[ix+1][iz].live==shotcells->live &&
				shotcells->cell[ix][iz+1].live==shotcells->live &&
				shotcells->cell[ix+1][iz+1].live==shotcells->live) {
				/* accumulate beam in cell */
				cellBeam(cells,shotcells,ix,iz);
			}
		}	
	}
}

static void cellBeam (Cells *cells, ShotCells *shotcells, int jx, int jz)
/*****************************************************************************
Accumulate Gaussian beam for one cell.
******************************************************************************
Input:
cells		pointer to cells
shotcells	pointer to shotcells
jx			x index of the cell in which to accumulate beam
jz			z index of the cell in which to accumulate beam
*****************************************************************************/
{
	int lx=cells->lx,lz=cells->lz,snx=cells->snx,nx=cells->nx,nz=cells->nz;
	float fx=cells->fx,dx=cells->dx;
	int fxn=fx/dx;
	float **g=cells->g;
	Cell **cell=cells->cell;
	Cell **cell1=shotcells->cell;
	BeamData *bd=cells->bd;
	int ntr=bd->ntr,nti=bd->nti;
	float dtr=bd->dtr,ftr=bd->ftr,dti=bd->dti,fti=bd->fti;
	complex **cf=bd->cf;
	int kxlo,kxhi,kzlo,kzhi,kx,kz,itr,iti;
	float odtr,odti,t00r,t01r,t10r,t11r,t00i,t01i,t10i,t11i,
		a00r,a01r,a10r,a11r,a00i,a01i,a10i,a11i,
		tx0r,tx1r,tx0i,tx1i,ax0r,ax1r,ax0i,ax1i,
		txzr,txzi,axzr,axzi,
		dtx0r,dtx0i,dtx1r,dtx1i,dax0r,dax0i,dax1r,dax1i,
		dtxzr,dtxzi,daxzr,daxzi,xdelta,zdelta,
		trn,tin,trfrac,mtrfrac,tifrac,mtifrac,
		cf0r,cf0i,cf1r,cf1i,cfr,cfi;
	complex *cf0,*cf1;

	/* inverse of time sampling intervals */
	odtr = 1.0/dtr;
	odti = 1.0/dti;

	/* complex time and amplitude for each corner */
	t00r = cell[jx][jz].tr+cell1[jx][jz].tr;
	t01r = cell[jx][jz+1].tr+cell1[jx][jz+1].tr;
	t10r = cell[jx+1][jz].tr+cell1[jx+1][jz].tr;
	t11r = cell[jx+1][jz+1].tr+cell1[jx+1][jz+1].tr;
	t00i = cell[jx][jz].ti+cell1[jx][jz].ti;
	t01i = cell[jx][jz+1].ti+cell1[jx][jz+1].ti;
	t10i = cell[jx+1][jz].ti+cell1[jx+1][jz].ti;
	t11i = cell[jx+1][jz+1].ti+cell1[jx+1][jz+1].ti;

	a00r = (cell[jx][jz].ar+cell1[jx][jz].ar)/2;
	a01r = (cell[jx][jz+1].ar+cell1[jx][jz+1].ar)/2;
	a10r = (cell[jx+1][jz].ar+cell1[jx+1][jz].ar)/2;
	a11r = (cell[jx+1][jz+1].ar+cell1[jx+1][jz+1].ar)/2;
	a00i = (cell[jx][jz].ai+cell1[jx][jz].ai)/2;
	a01i = (cell[jx][jz+1].ai+cell1[jx][jz+1].ai)/2;
	a10i = (cell[jx+1][jz].ai+cell1[jx+1][jz].ai)/2;
	a11i = (cell[jx+1][jz+1].ai+cell1[jx+1][jz+1].ai)/2;

	/* x and z samples for cell */
	kxlo = jx*lx+fxn;
	kxhi = kxlo+lx;
	if (kxhi>nx) kxhi =nx;
	kzlo = jz*lz;
	kzhi = kzlo+lz;
	if (kzhi>nz) kzhi = nz;

	/* fractional increments for linear interpolation */
	xdelta = 1.0/lx;
	zdelta = 1.0/lz;

	/* increments for times and amplitudes at top and bottom of cell */
	dtx0r = (t10r-t00r)*xdelta;
	dtx1r = (t11r-t01r)*xdelta;
	dtx0i = (t10i-t00i)*xdelta;
	dtx1i = (t11i-t01i)*xdelta;
	dax0r = (a10r-a00r)*xdelta;
	dax1r = (a11r-a01r)*xdelta;
	dax0i = (a10i-a00i)*xdelta;
	dax1i = (a11i-a01i)*xdelta;

	/* times and amplitudes at top-left and bottom-left of cell */
	tx0r = t00r;
	tx1r = t01r;
	tx0i = t00i;
	tx1i = t01i;
	ax0r = a00r;
	ax1r = a01r;
	ax0i = a00i;
	ax1i = a01i;

	/* loop over x samples */
	for (kx=kxlo; kx<kxhi; ++kx) {

		/* increments for time and amplitude */
		dtxzr = (tx1r-tx0r)*zdelta;
		dtxzi = (tx1i-tx0i)*zdelta;
		daxzr = (ax1r-ax0r)*zdelta;
		daxzi = (ax1i-ax0i)*zdelta;

		/* time and amplitude at top of cell */
		txzr = tx0r;
		txzi = tx0i;
		axzr = ax0r;
		axzi = ax0i;

		/* loop over z samples */
		for (kz=kzlo; kz<kzhi; ++kz) {

			/* index of imaginary time */
			iti = tin = (txzi-fti)*odti;
			if (iti<0 || iti>=nti-1) continue;

			/* pointers to left and right imaginary time samples */
			cf0 = cf[iti];
			cf1 = cf[iti+1];

			/* imaginary time linear interpolation coefficients */
			tifrac = tin-iti;
			mtifrac = 1.0-tifrac;

			/* index of real time */
			itr = trn = (txzr-ftr)*odtr;
			if (itr<0 || itr>=ntr-1) continue;

			/* real time linear interpolation coefficients */
			trfrac = trn-itr;
			mtrfrac = 1.0-trfrac;

			/* real and imaginary parts of complex beam data */
			cf0r = mtrfrac*cf0[itr].r+trfrac*cf0[itr+1].r;
			cf1r = mtrfrac*cf1[itr].r+trfrac*cf1[itr+1].r;
			cfr = mtifrac*cf0r+tifrac*cf1r;
			cf0i = mtrfrac*cf0[itr].i+trfrac*cf0[itr+1].i;
			cf1i = mtrfrac*cf1[itr].i+trfrac*cf1[itr+1].i;
			cfi = mtifrac*cf0i+tifrac*cf1i;

			/* accumulate beam */
			g[kx][kz] += axzr*cfr-axzi*cfi;

			/* increment time and amplitude */
			txzr += dtxzr;
			txzi += dtxzi;
			axzr += daxzr;
			axzi += daxzi;
		}

		/* increment times and amplitudes at top and bottom of cell */
		tx0r += dtx0r;
		tx1r += dtx1r;
		tx0i += dtx0i;
		tx1i += dtx1i;
		ax0r += dax0r;
		ax1r += dax1r;
		ax0i += dax0i;
	    ax1i += dax1i;
	}
}
