/*
 * upwindorderedMATHEVALv4
 * added in Ben Nolting's two boundary cases of bounce and default
 * replaced GNU libmatheval with expression_parser
 * 
 * system("R CMD SHLIB -I/usr/local/include -L/usr/local/lib expression_parser.c -lm")
 * system("R CMD SHLIB -I/usr/local/include -L/usr/local/lib upwindorderedv4.c expression_parser.o -lm")
 * currentupwindordered = "upwindorderedv4"
 * try( dyn.load(paste(currentupwindordered, ".so", sep = "")) )
 * 
 * upwindorderedMATHEVALv3
 * changed from static buffer size to dynamic buffer size for the equation strings
 * 
 *  upwindorderedMATHEVALv2
 * fix the assumption that LX1 and LY1 have to be negative

 * upwindorderedMATHEVAL
 * works with MATHEVAL library from GNU
 * can call the quasipotentialSTRINGSwithstoragePOINTERSreturnvalue function from within R
*/

/* Ordered Upwind Method */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include "expression_parser.h" 
#include <math.h>
#include <time.h>
#include <R.h>

#define ourPI 3.141592653589793
#define mabs(a) ((a) >= 0 ? (a) : -(a))
#define sgn(a) ((a) == 0 ? 0 : ((a) > 0  ? 1 : -1 ))
#define max(a,b) ((a) >= (b) ? (a) : (b))
#define min(a,b) ((a) <= (b) ? (a) : (b))
#define INFTY 1.0e+6 /* This is the "maximum" value of the quasi-potential. It is also the default value... everything is initialized to this values before the solution computation starts. So points with unknown values have value 1.0e+6 */
#define TOL 1.0e-12
#define BETA 0.0
#define BUFFERSIZE 1024	/* maximum size of the equation to be used */


/*Pull all of these out of define and make them global variables that 
 * get values given to C from R attached to them
 */

int NX = 260; /* This is the number of discrete units in the horizontal direction.*/
int NY = 260; /* This is the number of discrete units in the vertical direction. So the total number of discrete points on the grid is NX*NY*/

double LX1 = 5.0; /* lower bound v2 :::: NO LONGER APPLIES v2: The left bound of the domain is -LX1. So if you want to go from -5 to 10, LX1=5. Note that it is a good idea to give yourself space, and not set this equal to zero. The computation stops when you get to a boundary. */
double LX2 = 60.0; /* The right bound of the domain is LX2. So if you want to go from -5 to 10, LX2=10. Note that it is a good idea to give yourself space, and not set this equal to zero. The computation stops when you get to a boundary. */
double LY1 = 5.0; /* Lower bound v2 :::: NO LONGER APPLIES v2:Same idea as LX1, but for the vertical direction */
double LY2 = 60.0; /* Same idea as LX2, but for the vertical direction */

double FP1 = 6.60341; /* This is the x-coordinate of the initial point for the calculation (equilibirium or point on a limit cycle). */
double FP2 = 3.04537; /* This is the y-coordinate of the initial point for the calculation (equilibirium or point on a limit cycle). */

/* variables used by matheval */
char xbuffer[BUFFERSIZE];		/*stores the equation as a string*/
char ybuffer[BUFFERSIZE];
char *xbuff;
char *ybuff;
char *name[] = {"x", "y"};		/*variables used in the equations*/
char *filenamebuff;
double values[] = {0.0, 0.0};	/*values of variables used in equations*/

double bounceedge = 0.01;	/* used to compute a bounce edge for the bounce reflecting boundaries */


struct myvector {
    double x;
    double y;
} myvector;


struct mymatrix {
    double a11;
    double a12;
    double a21;
    double a22;
} mymatrix;

struct mysol {
    double g;
    char c;
} mysol;

int main(int argc, char **argv);
void quasipotential(double *storage, double *tempxmin, double *tempxmax, int *tempxsteps, double *tempymin, double *tempymax, int *tempysteps, double *tempeqx, double *tempeqy, char **equationx, int *lenequationx, char **equationy, int *lenequationy, char **tempfilename, int *templengthfilename, int *tempdatasave, char **tempchfield, double *tempbounceedge);


/* struct myvector myfield(double x,double y); /* B */
struct myvector myfieldchris(double x,double y);
double angle(double x,double y);
double length(double x,double y);
void param(void);
void ipoint(void);
void ordered_upwind(void);
struct mysol triangle_update(long ind,long ind0,long ind1);
double one_pt_update(long ind,long ind0);
void addtree(long ind); /* adds a node to the binary tree
                         of the "considered" points */
void updatetree(long ind); /* updates the binary tree */
void deltree(void); /* deletes the root of the binary tree */
struct mymatrix matrix_inverse(struct mymatrix matr);
struct mymatrix matrix_product(struct mymatrix a,struct mymatrix b);
struct mymatrix transpose(struct mymatrix matr);
struct myvector matr_vec(struct mymatrix matr,struct myvector vec);
double dotproduct(struct myvector a,struct myvector b);
struct myvector ga_plus_b(double lam,struct myvector a,struct myvector b);
double solve_quadratic(double a,double b,double c);
struct myvector getpoint(long ind);
/***************************************/

/* const long nx1=NX-1, ny1=NY-1, nxy=NX*NY; */
long nx1,ny1,nxy;
long count=0; /* # of considered points */
double h,hx,hy;
double *aB = 0; /* potential on the regular mesh */
struct myvector *B = 0; /* exp(2 beta v) */
int *ms = 0; /* Label for points in algorithm.  0 = 'Unknown', 1 = 'Considered', 2 = "in Accepted Front", 3 = "Accepted" */
double *g = 0; /* g is an array that stores the answers */
struct myvector *rcurr = 0; /* M grad g */
long *pos = 0; /* pos(index of mesh pt) = position in binary tree */
long *tree = 0; /* tree(position in the tree) = index of mesh pt */
double UPS, HUPS; /* unisotropy ratio. This measures how funky your vector field is. */
long *acf = 0; /* vector of indexes of points in the accepted front */
long *pacf = 0; /* position in the accepted front */
/* const long neii[8]={1, NX+1, NX, NX-1, -1, -NX-1, -NX, -NX+1 }; /* neighbor's indices */
long neii[8]; /*={1, NX+1, NX, NX-1, -1, -NX-1, -NX, -NX+1 }; /* neighbor's indices */
double xa, ya, xb, yb; /* potential minima */
char chfield = 'd';

/**************************************/
/**************************************/
/* *************************************
struct myvector myfieldchris(double x,double y) {
    struct myvector v;
    v.x = evaluator_evaluate_x_y(fx, x, y);
    v.y = evaluator_evaluate_x_y(fy, x, y);
    return v;
}
* **************************************/

/* variable_back used to read in strings
 * converts x and y to actual value in c code
*/
double chrisx;
double chrisy;

int variable_callback( void *user_data, const char *name, double *value ){
	// look up the variables by name
	if( strcmp( name, "x" ) == 0 ){
		// set return value, return true
		*value = chrisx;
		return PARSER_TRUE;
	} else if( strcmp( name, "y" ) == 0 ){
		// set return value, return true
		*value = chrisy;
		return PARSER_TRUE;
	}
	// failed to find variable, return false
	return PARSER_FALSE;
}

/* taken right from expression_parser::example.c
 * because function_callback does not automatically go to NULL
*/
int function_callback( void *user_data, const char *name, const int num_args, const double *args, double *value ){
	int i, max_args;
	double tmp;
	
	// example to show the user-data parameter, sets the maximum number of
	// arguments allowed for the following functions from the user-data function
	max_args = *((int*)user_data);
	
	if( strcmp( name, "max_value") == 0 && num_args >= 2 && num_args <= max_args ){
		// example 'maximum' function, returns the largest of the arguments, this and 
		// the min_value function implementation below allow arbitrary number of arguments
		tmp = args[0];
		for( i=1; i<num_args; i++ ){
			tmp = args[i] >= tmp ? args[i] : tmp;
		}
		// set return value and return true
		*value = tmp;
		return PARSER_TRUE;
	} else if( strcmp( name, "min_value" ) == 0 && num_args >= 2 && num_args <= max_args ){
		// example 'minimum' function, returns the smallest of the arguments
		tmp = args[0];
		for( i=1; i<num_args; i++ ){
			tmp = args[i] <= tmp ? args[i] : tmp;
		}
		// set return value and return true
		*value = tmp;
		return PARSER_TRUE;
	} 
	
	// failed to evaluate function, return false
	return PARSER_FALSE;
}

struct myvector myfieldchris(double x,double y) {
    struct myvector v;
/* hack until I can check for sure what variables I am using */ 
    chrisx = x;
    chrisy = y;
    int num_arguments = 4; /*need for parse_expression_with_callbacks; number picked and means nothing */
    switch( chfield ) {
        case 'p': /* positivevalues: Case to use if you want only positive values */
            if(x<0 && y<0)
            {v.x=1000000.0;
                v.y=1000000.0;}
            else if(x<0)
            {v.x=1000000.0;
                v.y=0.0;}
            else if(y<0)
            {v.x=0.0;
                v.y=100000.0;}
            else {
            	v.x = parse_expression_with_callbacks( xbuff, variable_callback, function_callback, &num_arguments );
                v.y = parse_expression_with_callbacks( ybuff, variable_callback, function_callback, &num_arguments );
            }
            break;
		case 'b': /* bounce : Case to use if you want reflecting boundaries */
            if(x<(LX1+bounceedge*NX*hx) && y<(LY1+bounceedge*NY*hy) )
            {v.x=1000000.0;
                v.y=1000000.0;}
            else if(x<(LX1+bounceedge*NX*hx))
            {v.x=1000000.0;
                v.y=0.0;}
            else if(y<(LY1+bounceedge*NY*hy))
            {v.x=0.0;
                v.y=100000.0;}
            else {
            	v.x = parse_expression_with_callbacks( xbuff, variable_callback, function_callback, &num_arguments );
                v.y = parse_expression_with_callbacks( ybuff, variable_callback, function_callback, &num_arguments );
            }
            break;
        case 'd': /*default */
			v.x = parse_expression_with_callbacks( xbuff, variable_callback, function_callback, &num_arguments );
			v.y = parse_expression_with_callbacks( ybuff, variable_callback, function_callback, &num_arguments );
            break;
        default:
			Rprintf("Should not be able to reach this in myfieldchris function");
            break;
/*            Rprintf("chfield = %c, please correct\n",chfield);
            exit(1);
            break;
*/
    }
    return v;
}

/*************************************/

void param() {
    long i,j,ind;
    double x,y;
    
    xa=FP1; ya=FP2;
/* v2: LX1, LY1 used to compute step size, 
 * so need to take range from L2 to L1 and divide by number of steps 
 * */
	hx=(LX2 - LX1)/(NX-1);
	Rprintf("hx = %g\n", hx);
	hy=(LY2-LY1)/(NY-1);
	Rprintf("hy = %g\n", hy);
    h=sqrt(hx*hx+hy*hy);
    for( i=0; i<NX; i++ ) {
/* v2: LX1 and LY1 used as beginning value, to compute "current" x and y
 * need to add stepsize times step number and add to lower bound
 * */
        x=LX1+hx*i;
        for( j=0; j<NY; j++ ) {
            ind=i+NX*j;
            y=LY1+hy*j;
            rcurr[ind].x=0.0;
            rcurr[ind].y=0.0;
            B[ind]=myfieldchris(x,y);
            aB[ind]=length(B[ind].x,B[ind].y);
            ms[ind]=0;
            g[ind]=INFTY;
        }
    }
}


/************************************/

void ipoint() {
    long i,j,ind,ind0,m,n,*indac;
	long *ic,nac=0,nc=0;;
    double gtemp;
    const long isur[4]={0, 1, NX+1, NX};
    struct myvector x,l;
    
/*    Rprintf("DECLARED PARAMETERS IN ipoint()\n"); */
    indac=(long *)malloc(4*sizeof(long));
/*    Rprintf("DECLARED indac IN ipoint()\n");*/

 /* v2: appears to convert the "current" value of x and y into an index value
  * so that (0+lower bound)/stepsize = 0 and (upper bound + lower bound)/ha = total number of steps
  * "current" value is lowest point, i.e. equilibrium 
  * (xa initialized to FP which is starting point, i.e. equilibrium)
  * */
    i=floor((xa-(LX1))/hx);
    j=floor((ya-(LY1))/hy);
/*    Rprintf("DECLARED i and j IN ipoint()\n");*/
    ind0=i+j*NX;
    for( m=0; m<4; m++ ) {
/*		Rprintf("DECLARED m = %li\n", m);*/
        ind=ind0+isur[m];
/*        Rprintf("DECLARED ind = %li\n", ind);*/
/*        Rprintf("BEFORE getpoint() CALL\n");*/
        x=getpoint(ind);
/*        Rprintf("AFTER getpoint() CALL\n");*/
        l.x=x.x-xa;
        l.y=x.y-ya;
/*        Rprintf("BEFORE dotproduct() CALL\n");*/
        gtemp=aB[ind]*length(l.x,l.y)-dotproduct(B[ind],l);
/*        Rprintf("AFTER dotproduct() CALL\n");*/
        g[ind]=min(g[ind],gtemp);
        if( ms[ind] == 0 ) {
            ms[ind]=1;
            addtree(ind);
        }
    }
}

/**********************************************/
/*** ordered upwind method ***/

void ordered_upwind(void) {
    long i,j,k,m,ind,ind0,ind1m,ind1p,ii,jj,ind1,indupdate,kx,ky,i0,j0,i1,j1,mycount=0;
    int nc; /* # of points that become considered at the current cycle */
    double x0,x1,y0,y1,g0,g1,x,y,len0,len1;
    long newind[8]; /* the indices of the new considered points */
    struct mysol sol;
    char newch, update;
    const long KX=20, KY=20;
    double gamma, bdotvec, det, avec, aux,a0,b0,a1,b1,HUPS;
    struct myvector vec, b, c, v0, v1;
    double xnewac, ynewac; /* x and y of the newly accepted point  */
    
	HUPS=max(KX,KY)*h;
    while( count > 0 ) {
        ind=tree[1];
        j=ind/NX;
        i=ind%NX;
/* v2: converts from index value to a float
 * lower bound + index value * step size 
 * */
        xnewac=LX1 + i*hx;
        ynewac=LY1 + j*hy;
        ms[ind]=2;
        deltree();
        mycount++;
        if( i==2 || i==nx1-2 || j==2 || j== ny1-2 || g[ind] >= INFTY-1) {
            Rprintf("%ld\t(%ld\t%ld) is accepted, g=%.4f\n",mycount,i,j,g[ind]);
            break; /* quit if we reach the boundary of the computational domain */
        }
        /*	Rprintf("%ld\t(%ld\t%ld) is accepted, g=%.4f\n",mycount,i,j,g[ind]); */
        
        /* update considered neighbors of the accepted point */
        for( k=0; k<8; k++ ) {
            ind1=ind+neii[k];
            if( ms[ind1] == 2 ) {
                ind1m=ind+neii[(k-1+8)%8];
                ind1p=ind+neii[(k+1)%8];
                if( ms[ind1m] >= 2 && ms[ind1p] >= 2 ) {
                    ms[ind1]=3;
                    /*	  Rprintf("ms[%ld, %ld] = 3\n",ind1%NX,ind1/NX); */
                }
                else {
                    g0=g[ind];
                    ii=ind1%NX;
                    jj=ind1/NX;
/* v2: converts index value to x or y value */
                    x1=LX1 + ii*hx;
                    y1=LY1 + jj*hy;
                    g1=g[ind1];
                    for( i0=max(i-KX,0); i0<=min(i+KX,nx1); i0++) for( j0=max(j-KY,0); j0<=min(j+KY,ny1); j0++ ) {
                        indupdate=i0+NX*j0;
                        update='y';
                        if( ms[indupdate] == 1 ) {
/* v2: converts index value to x or y value */
                            x=LX1 + i0*hx;
                            y=LY1 + j0*hy;
                            vec.x=x-xnewac;
                            vec.y=y-ynewac;
                            len0=length(vec.x,vec.y);
                            len1=length(x-x1,y-y1);
                            b=B[indupdate];
                            bdotvec=dotproduct(b,vec);
                            
                            if( update == 'y' && fabs(max(len0,len1)-length(xnewac-x1,ynewac-y1)-min(len0,len1)) > TOL ) {
                                sol=triangle_update(indupdate,ind,ind1);
                                g[indupdate]=min(g[indupdate],sol.g);
                                if( sol.g <= g[indupdate] ) {
                                    updatetree(indupdate);
                                    /*	    Rprintf("( %ld %ld ) has been updated, g = %.4f\n", i0,j0,g[indupdate]); */
                                }
                            }
                        }
                    }
                }
            }
        }
        
        nc=0;
        for( k=0;k<8;k++) {
            ind1=ind+neii[k];
            if( ms[ind1]==0 )  {
                newind[nc]=ind1;
                nc++;
            }
        }
        for( m=0; m<nc; m++ ) {
            indupdate=newind[m];
            newch='n';
            i=indupdate%NX;
            j=indupdate/NX;
            
/* v2: converts index values to x or y values */
            x=LX1 + i*hx;
            y=LY1 + j*hy;
            vec.x=x-xnewac;
            vec.y=y-ynewac;
            b=B[indupdate];
            bdotvec=dotproduct(b,vec);
            
            /*   Rprintf("New Considered point: (%ld %ld)\n",i,j); */
            for( i0=max(0,i-KX); i0<=min(nx1,i+KX); i0++ ) for( j0=max(j-KY,0); j0<=min(ny1,j+KY); j0++ ) {
                ind0=i0+NX*j0;
/* v2: converts index values to x or y values */
                x0=LX1 + hx*i0;
                y0=LY1 + hy*j0;
                if( ms[ind0] == 2 || (ms[ind0]==3 && fabs(i-i0)<1.5 && fabs(j-j0)<1.5) ) {
                    update='y';
                    if( bdotvec > 0.0 ) {
                        v0.x=x-x0;
                        v0.y=y-y0;
                        a0=(v0.x*b.y-v0.y*b.x)/det;
                        b0=(v0.y*c.x-v0.x*c.y)/det;
                        if( a0 <= 0.0 || b0 <= 0.0 ) update='n';
                    }
                    for( k=0; k<8; k++ ) {
                        ind1=ind0+neii[k];
                        if( ms[ind1] == 2 ) {
                            update='y';
                            i1=ind1%NX;
                            j1=ind1/NX;
/* v2: converts index values to x or y values */
                            x1=LX1 + hx*i1;
                            y1=LY1 + hy*j1;
                            if( update == 'n' ) {
                                v1.x=x-x1;
                                v1.y=y-y1;
                                a1=(v1.x*b.y-v1.y*b.x)/det;
                                b1=(v1.y*c.x-v1.x*c.y)/det;
                                if( a1 >= 0.0 && b1 >= 0.0 ) update ='y';
                            }
                            len0=length(x-x0,y-y0);
                            len1=length(x-x1,y-y1);
                            if( bdotvec < 0.0 ) {
                                if( min(len0,len1) > HUPS ) update = 'n';
                            }
                            if( update == 'y' &&  fabs(max(len0,len1)-length(x0-x1,y0-y1)-min(len0,len1)) > TOL ) {
                                sol=triangle_update(indupdate,ind0,ind1);
                                g[indupdate]=min(g[indupdate],sol.g);
                                newch='y';
                            }
                        }
                    }
                }
            }
            if( newch == 'y' ) {
                addtree(indupdate);
                ms[indupdate]=1;
            }
            /*  Rprintf("   ( %ld\t%ld ) becomes Considered, g=%.4f\n",i,j,g[indupdate]); */
        }
        
    } /* end while ( count > 0 ) */
}




/*********************************************/
/*** triangle update ***/

struct mysol triangle_update(long ind,long ind0,long ind1) {
    struct myvector x,x0,x1,l0,l1;
    struct myvector p0,p1,a,b,qa,qb,coefs,vec,gradg,phis;
    struct mymatrix pmatr,pinv,pinvT;
    struct mymatrix qmatr;
    double d0,d1,aeq,beq,ceq,myf;
    double gtent,gtent0,gtent1;
    struct mysol sol;
    
    x=getpoint(ind);
    x0=getpoint(ind0);
    x1=getpoint(ind1);
    l0.x=x.x-x0.x;
    l0.y=x.y-x0.y;
    l1.x=x.x-x1.x;
    l1.y=x.y-x1.y;
    d0=length(l0.x,l0.y);
    d1=length(l1.x,l1.y);
    p0.x=l0.x/d0;
    p0.y=l0.y/d0;
    p1.x=l1.x/d1;
    p1.y=l1.y/d1;
    pmatr.a11=p0.x;
    pmatr.a12=p0.y;
    pmatr.a21=p1.x;
    pmatr.a22=p1.y;
    pinv=matrix_inverse(pmatr);
    qmatr=matrix_inverse(matrix_product(pmatr,transpose(pmatr)));
    a.x=1/d0;
    a.y=1/d1;
    b.x=-g[ind0]/d0;
    b.y=-g[ind1]/d1;
    qa=matr_vec(qmatr,a);
    qb=matr_vec(qmatr,b);
    /* compute coefficients for the quadratic equation for g(x,y) */
    aeq=dotproduct(a,qa);
    beq=2.0*(dotproduct(b,qa)+dotproduct(B[ind],matr_vec(pinv,a)));
    ceq=dotproduct(b,qb)+2.0*dotproduct(matr_vec(pinv,b),B[ind]);
    gtent=solve_quadratic(aeq,beq,ceq);
    if( gtent < INFTY ) {
        /* check upwinding criterium */
        /* compute coefs=P^{-T}(P^{-1}(gtent*a+b)+B)), check coefs.x>0, coefs.y>0 */
        vec=ga_plus_b(gtent,a,b);
        pinvT=transpose(pinv);
        gradg=matr_vec(pinv,vec);
        phis.x=(gradg.x+B[ind].x)/aB[ind];
        phis.y=(gradg.y+B[ind].y)/aB[ind];
        coefs=matr_vec(pinvT,phis);
        if( coefs.x >= 0.0 && coefs.y >= 0.0 ) {
            sol.g=gtent;
            sol.c='y';
        }
        else {
            sol.g=INFTY;
            sol.c='n';
        }
    }
    else sol.c='n';
        if( sol.c == 'n' ) {
            gtent0=one_pt_update(ind,ind0);
            gtent1=one_pt_update(ind,ind1);
            sol.g=min(gtent0,gtent1);
            sol.c='y';
        }
    return sol;
}

/*-----------*/

double one_pt_update(long ind,long ind0) {
    struct myvector x,x0,l;
    double gtemp;
    
    x=getpoint(ind);
    x0=getpoint(ind0);
    l.x=x.x-x0.x;
    l.y=x.y-x0.y;
    gtemp=g[ind0]+aB[ind]*length(l.x,l.y)-dotproduct(B[ind],l);
    return gtemp;
}

/*-------------*/

struct myvector getpoint(long ind) {
    struct myvector l;
/*	Rprintf("ALL VALUES IN GETPOINT: %li %g %i %g %g %g\n", ind, hx, NX, LX1, hy, LY1); */
/* v2: takes an index value and converts it to an x,y value
 * The reason this looks strange is that the storage matrix g is a 1 by NX*NY array
 * as opposed to a two-dimensional array
 * therefore the % and / and the use of NX in both cases.
 * */
    l.x = LX1 + hx*(ind%NX) ;
    l.y = LY1 + hy*(ind/NX);
    return l;
}

/**********************************************/
/*** linear algebra ***/

double angle(double x,double y) {
    double ang;
    if( y >= 0.0 ) ang=acos(x/sqrt(x*x+y*y));
    else ang=2.0*ourPI-acos(x/sqrt(x*x+y*y));
    return ang;
}

double length(double x,double y) {
    return sqrt(x*x+y*y);
}

struct mymatrix matrix_inverse(struct mymatrix matr) {
    struct mymatrix mi;
    double rdet;
    
    rdet=1.0/(matr.a11*matr.a22-matr.a12*matr.a21);
    mi.a11=matr.a22*rdet;
    mi.a12=-matr.a12*rdet;
    mi.a21=-matr.a21*rdet;
    mi.a22=matr.a11*rdet;
    return mi;
}

struct mymatrix matrix_product(struct mymatrix a,struct mymatrix b) {
    struct mymatrix c;
    
    c.a11=a.a11*b.a11+a.a12*b.a21;
    c.a12=a.a11*b.a12+a.a12*b.a22;
    c.a21=a.a21*b.a11+a.a22*b.a21;
    c.a22=a.a21*b.a12+a.a22*b.a22;
    return c;
}

struct mymatrix transpose(struct mymatrix matr) {
    struct mymatrix mt;
    
    mt.a11=matr.a11;
    mt.a22=matr.a22;
    mt.a12=matr.a21;
    mt.a21=matr.a12;
    return mt;
}

struct myvector matr_vec(struct mymatrix matr,struct myvector vec) {
    struct myvector c;
    
    c.x=matr.a11*vec.x+matr.a12*vec.y;
    c.y=matr.a21*vec.x+matr.a22*vec.y;
    return c;
}

double dotproduct(struct myvector a,struct myvector b) {
    return a.x*b.x+a.y*b.y;
}

struct myvector ga_plus_b(double lam,struct myvector a,struct myvector b) {
    struct myvector c;
    
    c.x=lam*a.x+b.x;
    c.y=lam*a.y+b.y;
    return c;
}

double solve_quadratic(double a,double b,double c) {
    double discr,sol;
    
    discr=b*b-4.0*a*c;
    if( discr < 0.0 ) sol=INFTY;
    else sol=0.5*(-b+sqrt(discr))/a;
    return sol;
}

/**************************************************************/
/************ FUNCTIONS RELATED TO THE BINARY TREE ***************/

void addtree(long ind) {
    long loc, ptemp;
    long indp, indc;
    char ch;
    
    count++;
    tree[count]=ind;
    pos[ind]=count;
    if( count > 1 ) {
        loc=count;
        indc=tree[loc];
        indp=tree[loc/2];
        ch=( g[indc] < g[indp] ) ? 'y' : 'n';
        while( ch == 'y' ) {
            ptemp=pos[indc];
            pos[indc]=pos[indp];
            tree[loc/2]=indc;
            pos[indp]=ptemp;
            tree[loc]=indp;
            loc=loc/2;
            if( loc > 1 ) {
                indc=tree[loc];
                indp=tree[loc/2];
                ch=( g[indc] < g[indp] ) ? 'y' : 'n';
            }
            else ch='n';
        }
    }
}

/*------------------------------------------------------------------*/

void updatetree(long ind) {
    long loc, lcc;
    double g0,g1,g2;
    
    g0=g[ind];
    loc=pos[ind];
    while( loc > 1 && g0 < g[tree[loc/2]] ) {
        tree[loc]=tree[loc/2];
        pos[tree[loc]]=loc;
        loc=loc/2;
        tree[loc]=ind;
        pos[tree[loc]]=loc;
    }
    g1=g[tree[loc*2]];
    g2=g[tree[loc*2+1]];
    lcc=count;
    while( (loc*2 <= count && g0 > g[tree[loc*2]]) || (loc*2+1 <= count && g0 > g[tree[loc*2+1]]) )  {
        lcc=( loc*2+1 <=count && g[tree[loc*2+1]] < g[tree[loc*2]] ) ? loc*2+1 : loc*2;
        tree[loc]=tree[lcc];
        pos[tree[loc]]=loc;
        loc=lcc;
        tree[loc]=ind;
        pos[tree[loc]]=loc;
    }
}

/*---------------------------------------------------------------------*/


/* deletes root of the binary tree */
void deltree() {
    long loc, ptemp, ind, lcc, ic, ic1, ic2, mind;
    char chd, ch='n';;
    
    mind=tree[1];
    pos[tree[1]]=0;
    tree[1]=tree[count];
    pos[tree[1]]=1;
    count--;
    loc=1;
    ind=tree[1];
    lcc=2*loc;
    if( lcc < count )  {
        ic1=tree[lcc];
        ic2=tree[lcc+1];
        if( (g[ind]) > (g[ic1]) || (g[ind]) > (g[ic2]) ) {
            if( (g[ic1]) <= (g[ic2]) )  {
                chd='l';
                ic=ic1;
            }
            else {
                chd='r';
                ic=ic2;
                lcc++;
            }
        }
        else chd='n';
    }
    else if( lcc == count ) {
        ic=tree[lcc];
        if( (g[ind]) > (g[ic]) ) {chd='l'; if(ch=='y') Rprintf("left\n");}
        else chd='n';
    }
    else chd='n';
    while( chd != 'n' ) {
        ptemp=pos[ind];
        pos[ind]=pos[ic];
        tree[loc]=ic;
        pos[ic]=ptemp;
        tree[lcc]=ind;
        loc=lcc;
        lcc=2*loc;
        if( lcc < count )  {
            ic1=tree[lcc];
            ic2=tree[lcc+1];
            if( (g[ind]) > (g[ic1]) || (g[ind]) > (g[ic2]) ) {
                if( (g[ic1]) <= (g[ic2]) )  {
                    chd='l';
                    ic=ic1;
                }
                else {
                    chd='r';
                    ic=ic2;
                    lcc++;
                }
            }
            else chd='n';
        }
        else if( lcc == count ) {
            ic=tree[lcc];
            if(ch=='y') Rprintf("child: loc(%li)=%li, t1=%.12e\n",ic1,lcc,g[ic1]);
            if( (g[ind]) > (g[ic]) ) { chd='l';if(ch=='y') Rprintf("left\n");}
            else chd='n';
        }
        else chd='n';
    } /* end while( chd != 'n' ) */
}


/********************************************************/
void quasipotential(double *storage, double *tempxmin, double *tempxmax, int *tempxsteps, double *tempymin, double *tempymax, int *tempysteps, double *tempeqx, double *tempeqy, char **equationx, int *lenequationx, char **equationy, int *lenequationy, char **tempfilename, int *templengthfilename, int *tempdatasave, char **tempchfield, double *tempbounceedge) {
/* Assign function parameter values to variables defined in C code*/
/* x range, y range, and starting values */
	LX1 = tempxmin[0]; LX2 = tempxmax[0]; NX = tempxsteps[0]; 
	LY1 = tempymin[0]; LY2 = tempymax[0]; NY = tempysteps[0];
	FP1 = tempeqx[0]; 
	FP2 = tempeqy[0];
/* variables for edge handling */
	char *temptempchfield = *tempchfield;
	chfield = temptempchfield[0]; 
	bounceedge = tempbounceedge[0];
	switch( chfield ) {
        case 'p':
			break;
		case 'b':
			break;
		case 'd':
			break;
		default:
			Rprintf("chfield must be (d)efauly, (p)ositive, or (b)ounce\n");
			return;
	}
			

/* handling rhs of equations */
	int lengthofequationx = lenequationx[0];
	int lengthofequationy = lenequationy[0];
	xbuff = malloc(sizeof(char)*(lengthofequationx+1));
	ybuff = malloc(sizeof(char)*(lengthofequationy+1));
	
	char *tempchrisx = *equationx;
	char *tempchrisy = *equationy;
	for (int jchris = 0; jchris < lengthofequationx; jchris++) {
		// xbuffer[jchris] = **equationx[jchris];
		/*xbuffer[jchris] = tempchrisx[jchris];	/* for fixed size buffer*/
		xbuff[jchris] = tempchrisx[jchris];
	} 
/*	xbuffer[lengthofequationx] = '\0';	*/
	xbuff[lengthofequationx] = '\0';
	
	for (int jchris = 0; jchris < lengthofequationy; jchris++) {
		// ybuffer[jchris] = **equationy[jchris];
		/*ybuffer[jchris] = tempchrisy[jchris]; /* for fixed size buffer */
		ybuff[jchris] = tempchrisy[jchris];
	} 
	/*ybuffer[lengthofequationy] = '\0';	*/
	ybuff[lengthofequationy] = '\0';
	
/* Determine filename if file saved to harddrive */
	char *filename;
	int datasave = tempdatasave[0]; /*what type of save do we need to do? data to R, data to HD, data to both */
	Rprintf("Creating file name.\n");
	int lengthfilename = templengthfilename[0];
	if (lengthfilename == 0) {
	/* use default naming scheme */
		char tempfilename[200];
		filename = tempfilename;
		char tempname[10];
		strcpy(filename,"defaultname-");
		sprintf(tempname, "x%4.4f", FP1);
		strcat(filename, tempname);
		sprintf(tempname, "y%4.4f", FP2);
		strcat(filename, tempname);
		strcat(filename, ".txt");
	} else {
		filenamebuff = malloc(sizeof(char)*(lengthfilename+1));
		char *tempchrisfilename = *tempfilename;
		for (int jchris = 0; jchris < lengthfilename; jchris++) {
			filenamebuff[jchris] = tempchrisfilename[jchris];
		}
		filenamebuff[lengthfilename] = '\0';
		filename = filenamebuff;
	}
	Rprintf("File name created.\n");

/* The workhorse */
	nx1=NX-1; ny1=NY-1; nxy=NX*NY;
	neii[0]=1; neii[1]=NX+1; neii[2]=NX; neii[3]=NX-1; 
	neii[4]=-1; neii[5]=-NX-1; neii[6]=-NX; neii[7]=-NX+1;
	
	aB = malloc(sizeof(double)*NX*NY); 
	B = malloc(sizeof(myvector)*NX*NY); 
	ms = malloc(sizeof(int)*NX*NY); 
	g = malloc(sizeof(double)*NX*NY); 
	rcurr = malloc(sizeof(myvector)*NX*NY);
	pos = malloc(sizeof(long)*NX*NY); 
	tree = malloc(sizeof(long)*NX*NY);
	acf = malloc(sizeof(long)*2*(NX+NY)); 
	pacf = malloc(sizeof(long)*NX*NY);
	
	Rprintf("Completed Memory Allocation\n");
    long i,j,ind;
    clock_t CPUbegin;
    double cpu;
    
    Rprintf("equationx = %s\n", xbuff);
	Rprintf("equationy = %s\n", ybuff);
    
    param();
    Rprintf("Finished Loading Parameters\n");
    CPUbegin=clock();
    ipoint();
    Rprintf("Finished ipoint() function\n");
    ordered_upwind();
    cpu=(clock()-CPUbegin)/((double)CLOCKS_PER_SEC);
    Rprintf("Finished ordered_upwind() function\n");
    Rprintf("cputime = %g\n",cpu);
    
	FILE *fg;
/* Write data to use some where, some how */
	switch(datasave) {
	case 1: /* does not save to R, only saves to hard drive */
		fg=fopen(filename, "w");
		Rprintf("File opened.\n");
		Rprintf("In datasave case 1\n");
    
		ind=0;
		for( j=0; j<(NY); j++ ) {
			for( i=0; i<(NX-1); i++ ) {
				fprintf(fg,"%.4e\t",g[ind]);
				ind++;
			}
			fprintf(fg,"%.4e",g[(ind)]);
			ind++;
			fprintf(fg,"\n");
		}
		fclose(fg);
		break;
	case 2: /* saves to R, but does not save to hard drive */
		Rprintf("Saves only to R\n");
		Rprintf("In datasave case 2\n");
		ind=0;
		for( j=0; j<(NY); j++ ) {
			for( i=0; i<(NX-1); i++ ) {
				storage[ind] = g[ind];
				ind++;
			}
			storage[ind] = g[ind];
			ind++;
		}
		break;
	case 3:	/* saves to R and saves to hard drive */
		fg=fopen(filename, "w");
		Rprintf("In datasave case 3\n");
		Rprintf("File opened.\n");
    
		ind=0;
		for( j=0; j<(NY); j++ ) {
			for( i=0; i<(NX-1); i++ ) {
				fprintf(fg,"%.4e\t",g[ind]);
				storage[ind] = g[ind];
				ind++;
			}
			fprintf(fg,"%.4e",g[(ind)]);
			storage[ind] = g[ind];
			ind++;
			fprintf(fg,"\n");
		}
		fclose(fg);
		break;
	default:
		Rprintf("Running testrun.  You are not saving any data.\n");
		break;
	}
	
	


/*Free everything, although I think I missed stuff */
	if (lengthfilename != 0) {free(filenamebuff);}
	free(xbuff); free(ybuff); 
    free(aB); free(B); free(ms); free(g); free(rcurr); free(pos); free(tree); free(acf); free(pacf); 

	Rprintf("Successful.  Exiting C code\n");
}



/*** main ***/
int main (int argc, char **argv) {
	Rprintf("This code has been stripped of everything.  Only useful by R.");
	return 0;
}
