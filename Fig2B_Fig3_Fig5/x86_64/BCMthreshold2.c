/* Created by Language version: 7.7.0 */
/* VECTORIZED */
#define NRN_VECTORIZED 1
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "scoplib_ansi.h"
#undef PI
#define nil 0
#include "md1redef.h"
#include "section.h"
#include "nrniv_mf.h"
#include "md2redef.h"
 
#if METHOD3
extern int _method3;
#endif

#if !NRNGPU
#undef exp
#define exp hoc_Exp
extern double hoc_Exp(double);
#endif
 
#define nrn_init _nrn_init__BCMthreshold2
#define _nrn_initial _nrn_initial__BCMthreshold2
#define nrn_cur _nrn_cur__BCMthreshold2
#define _nrn_current _nrn_current__BCMthreshold2
#define nrn_jacob _nrn_jacob__BCMthreshold2
#define nrn_state _nrn_state__BCMthreshold2
#define _net_receive _net_receive__BCMthreshold2 
#define state state__BCMthreshold2 
 
#define _threadargscomma_ _p, _ppvar, _thread, _nt,
#define _threadargsprotocomma_ double* _p, Datum* _ppvar, Datum* _thread, _NrnThread* _nt,
#define _threadargs_ _p, _ppvar, _thread, _nt
#define _threadargsproto_ double* _p, Datum* _ppvar, Datum* _thread, _NrnThread* _nt
 	/*SUPPRESS 761*/
	/*SUPPRESS 762*/
	/*SUPPRESS 763*/
	/*SUPPRESS 765*/
	 extern double *getarg();
 /* Thread safe. No static _p or _ppvar. */
 
#define t _nt->_t
#define dt _nt->_dt
#define d0 _p[0]
#define p0 _p[1]
#define scounttau _p[2]
#define alpha _p[3]
#define scount0 _p[4]
#define bcm_scale _p[5]
#define alpha_scount _p[6]
#define d _p[7]
#define p _p[8]
#define tspike _p[9]
#define scount _p[10]
#define flagOLD _p[11]
#define boltzfactor _p[12]
#define pf _p[13]
#define output _p[14]
#define Dscount _p[15]
#define v _p[16]
#define _g _p[17]
#define _tsav _p[18]
#define _nd_area  *_ppvar[0]._pval
 
#if MAC
#if !defined(v)
#define v _mlhv
#endif
#if !defined(h)
#define h _mlhh
#endif
#endif
 
#if defined(__cplusplus)
extern "C" {
#endif
 static int hoc_nrnpointerindex =  -1;
 static Datum* _extcall_thread;
 static Prop* _extcall_prop;
 /* external NEURON variables */
 /* declaration of user functions */
 static int _mechtype;
extern void _nrn_cacheloop_reg(int, int);
extern void hoc_register_prop_size(int, int, int);
extern void hoc_register_limits(int, HocParmLimits*);
extern void hoc_register_units(int, HocParmUnits*);
extern void nrn_promote(Prop*, int, int);
extern Memb_func* memb_func;
 
#define NMODL_TEXT 1
#if NMODL_TEXT
static const char* nmodl_file_text;
static const char* nmodl_filename;
extern void hoc_reg_nmodl_text(int, const char*);
extern void hoc_reg_nmodl_filename(int, const char*);
#endif

 extern Prop* nrn_point_prop_;
 static int _pointtype;
 static void* _hoc_create_pnt(_ho) Object* _ho; { void* create_point_process();
 return create_point_process(_pointtype, _ho);
}
 static void _hoc_destroy_pnt();
 static double _hoc_loc_pnt(_vptr) void* _vptr; {double loc_point_process();
 return loc_point_process(_pointtype, _vptr);
}
 static double _hoc_has_loc(_vptr) void* _vptr; {double has_loc_point();
 return has_loc_point(_vptr);
}
 static double _hoc_get_loc_pnt(_vptr)void* _vptr; {
 double get_loc_point_process(); return (get_loc_point_process(_vptr));
}
 extern void _nrn_setdata_reg(int, void(*)(Prop*));
 static void _setdata(Prop* _prop) {
 _extcall_prop = _prop;
 }
 static void _hoc_setdata(void* _vptr) { Prop* _prop;
 _prop = ((Point_process*)_vptr)->_prop;
   _setdata(_prop);
 }
 /* connect user functions to hoc names */
 static VoidFunc hoc_intfunc[] = {
 0,0
};
 static Member_func _member_func[] = {
 "loc", _hoc_loc_pnt,
 "has_loc", _hoc_has_loc,
 "get_loc", _hoc_get_loc_pnt,
 0, 0
};
 /* declare global and static user variables */
#define boltzman boltzman_BCMthreshold2
 double boltzman = 0;
 /* some parameters have upper and lower limits */
 static HocParmLimits _hoc_parm_limits[] = {
 0,0,0
};
 static HocParmUnits _hoc_parm_units[] = {
 "tspike", "ms",
 0,0
};
 static double delta_t = 0.01;
 /* connect global user variables to hoc */
 static DoubScal hoc_scdoub[] = {
 "boltzman_BCMthreshold2", &boltzman_BCMthreshold2,
 0,0
};
 static DoubVec hoc_vdoub[] = {
 0,0,0
};
 static double _sav_indep;
 static void nrn_alloc(Prop*);
static void  nrn_init(_NrnThread*, _Memb_list*, int);
static void nrn_state(_NrnThread*, _Memb_list*, int);
 static void nrn_cur(_NrnThread*, _Memb_list*, int);
static void  nrn_jacob(_NrnThread*, _Memb_list*, int);
 
#define _watch_array _ppvar + 3 
 static void _hoc_destroy_pnt(_vptr) void* _vptr; {
   Prop* _prop = ((Point_process*)_vptr)->_prop;
   if (_prop) { _nrn_free_watch(_prop->dparam, 3, 2);}
   destroy_point_process(_vptr);
}
 
static int _ode_count(int);
static void _ode_map(int, double**, double**, double*, Datum*, double*, int);
static void _ode_spec(_NrnThread*, _Memb_list*, int);
static void _ode_matsol(_NrnThread*, _Memb_list*, int);
 
#define _cvode_ieq _ppvar[5]._i
 static void _ode_matsol_instance1(_threadargsproto_);
 /* connect range variables in _p that hoc is supposed to know about */
 static const char *_mechanism[] = {
 "7.7.0",
"BCMthreshold2",
 "d0",
 "p0",
 "scounttau",
 "alpha",
 "scount0",
 "bcm_scale",
 0,
 "alpha_scount",
 "d",
 "p",
 "tspike",
 0,
 "scount",
 0,
 0};
 
extern Prop* need_memb(Symbol*);

static void nrn_alloc(Prop* _prop) {
	Prop *prop_ion;
	double *_p; Datum *_ppvar;
  if (nrn_point_prop_) {
	_prop->_alloc_seq = nrn_point_prop_->_alloc_seq;
	_p = nrn_point_prop_->param;
	_ppvar = nrn_point_prop_->dparam;
 }else{
 	_p = nrn_prop_data_alloc(_mechtype, 19, _prop);
 	/*initialize range parameters*/
 	d0 = 0;
 	p0 = 0;
 	scounttau = 0;
 	alpha = 0;
 	scount0 = 0;
 	bcm_scale = 0.2;
  }
 	_prop->param = _p;
 	_prop->param_size = 19;
  if (!nrn_point_prop_) {
 	_ppvar = nrn_prop_datum_alloc(_mechtype, 6, _prop);
  }
 	_prop->dparam = _ppvar;
 	/*connect ionic variables to this model*/
 
}
 static void _initlists();
  /* some states have an absolute tolerance */
 static Symbol** _atollist;
 static HocStateTolerance _hoc_state_tol[] = {
 0,0
};
 
#define _tqitem &(_ppvar[2]._pvoid)
 static void _net_receive(Point_process*, double*, double);
 static void _net_init(Point_process*, double*, double);
 extern Symbol* hoc_lookup(const char*);
extern void _nrn_thread_reg(int, int, void(*)(Datum*));
extern void _nrn_thread_table_reg(int, void(*)(double*, Datum*, Datum*, _NrnThread*, int));
extern void hoc_register_tolerance(int, HocStateTolerance*, Symbol***);
extern void _cvode_abstol( Symbol**, double*, int);

 void _BCMthreshold2_reg() {
	int _vectorized = 1;
  _initlists();
 	_pointtype = point_register_mech(_mechanism,
	 nrn_alloc,nrn_cur, nrn_jacob, nrn_state, nrn_init,
	 hoc_nrnpointerindex, 1,
	 _hoc_create_pnt, _hoc_destroy_pnt, _member_func);
 _mechtype = nrn_get_mechtype(_mechanism[1]);
     _nrn_setdata_reg(_mechtype, _setdata);
 #if NMODL_TEXT
  hoc_reg_nmodl_text(_mechtype, nmodl_file_text);
  hoc_reg_nmodl_filename(_mechtype, nmodl_filename);
#endif
  hoc_register_prop_size(_mechtype, 19, 6);
  hoc_register_dparam_semantics(_mechtype, 0, "area");
  hoc_register_dparam_semantics(_mechtype, 1, "pntproc");
  hoc_register_dparam_semantics(_mechtype, 2, "netsend");
  hoc_register_dparam_semantics(_mechtype, 3, "watch");
  hoc_register_dparam_semantics(_mechtype, 4, "watch");
  hoc_register_dparam_semantics(_mechtype, 5, "cvodeieq");
 	hoc_register_cvode(_mechtype, _ode_count, _ode_map, _ode_spec, _ode_matsol);
 	hoc_register_tolerance(_mechtype, _hoc_state_tol, &_atollist);
 pnt_receive[_mechtype] = _net_receive;
 pnt_receive_init[_mechtype] = _net_init;
 pnt_receive_size[_mechtype] = 1;
 	hoc_register_var(hoc_scdoub, hoc_vdoub, hoc_intfunc);
 	ivoc_help("help ?1 BCMthreshold2 /Users/shkim/Research/neurcomp/Neural-information-processing-and-computations-of-two-input-synapses/Fig2BandFig3/x86_64/BCMthreshold2.mod\n");
 hoc_register_limits(_mechtype, _hoc_parm_limits);
 hoc_register_units(_mechtype, _hoc_parm_units);
 }
static int _reset;
static char *modelname = "";

static int error;
static int _ninits = 0;
static int _match_recurse=1;
static void _modl_cleanup(){ _match_recurse=1;}
 
static int _ode_spec1(_threadargsproto_);
/*static int _ode_matsol1(_threadargsproto_);*/
 static int _slist1[1], _dlist1[1];
 static int state(_threadargsproto_);
 
/*CVODE*/
 static int _ode_spec1 (double* _p, Datum* _ppvar, Datum* _thread, _NrnThread* _nt) {int _reset = 0; {
   Dscount = - scount / pf ;
   }
 return _reset;
}
 static int _ode_matsol1 (double* _p, Datum* _ppvar, Datum* _thread, _NrnThread* _nt) {
 Dscount = Dscount  / (1. - dt*( ( - 1.0 ) / pf )) ;
  return 0;
}
 /*END CVODE*/
 static int state (double* _p, Datum* _ppvar, Datum* _thread, _NrnThread* _nt) { {
    scount = scount + (1. - exp(dt*(( - 1.0 ) / pf)))*(- ( 0.0 ) / ( ( - 1.0 ) / pf ) - scount) ;
   }
  return 0;
}
 
static double _watch1_cond(_pnt) Point_process* _pnt; {
 	double* _p; Datum* _ppvar; Datum* _thread; _NrnThread* _nt;
	_thread= (Datum*)0; _nt = (_NrnThread*)_pnt->_vnt;
 	_p = _pnt->_prop->param; _ppvar = _pnt->_prop->dparam;
	v = NODEV(_pnt->node);
	return  ( v ) - ( 0.0 ) ;
}
 
static void _net_receive (_pnt, _args, _lflag) Point_process* _pnt; double* _args; double _lflag; 
{  double* _p; Datum* _ppvar; Datum* _thread; _NrnThread* _nt;
   int _watch_rm = 0;
   _thread = (Datum*)0; _nt = (_NrnThread*)_pnt->_vnt;   _p = _pnt->_prop->param; _ppvar = _pnt->_prop->dparam;
  if (_tsav > t){ extern char* hoc_object_name(); hoc_execerror(hoc_object_name(_pnt->ob), ":Event arrived out of order. Must call ParallelContext.set_maxstep AFTER assigning minimum NetCon.delay");}
 _tsav = t;   if (_lflag == 1. ) {*(_tqitem) = 0;}
 {
   if ( _lflag  == 0.0 ) {
     output = 0.0 ;
     }
   else if ( _lflag  == 2.0 ) {
     tspike = t ;
     output = 1.0 ;
     }
   else {
       _nrn_watch_activate(_watch_array, _watch1_cond, 1, _pnt, _watch_rm++, 2.0);
 }
   } }
 
static void _net_init(Point_process* _pnt, double* _args, double _lflag) {
       double* _p = _pnt->_prop->param;
    Datum* _ppvar = _pnt->_prop->dparam;
    Datum* _thread = (Datum*)0;
    _NrnThread* _nt = (_NrnThread*)_pnt->_vnt;
 _args[0] = 0.0 ;
   }
 
static int _ode_count(int _type){ return 1;}
 
static void _ode_spec(_NrnThread* _nt, _Memb_list* _ml, int _type) {
   double* _p; Datum* _ppvar; Datum* _thread;
   Node* _nd; double _v; int _iml, _cntml;
  _cntml = _ml->_nodecount;
  _thread = _ml->_thread;
  for (_iml = 0; _iml < _cntml; ++_iml) {
    _p = _ml->_data[_iml]; _ppvar = _ml->_pdata[_iml];
    _nd = _ml->_nodelist[_iml];
    v = NODEV(_nd);
     _ode_spec1 (_p, _ppvar, _thread, _nt);
 }}
 
static void _ode_map(int _ieq, double** _pv, double** _pvdot, double* _pp, Datum* _ppd, double* _atol, int _type) { 
	double* _p; Datum* _ppvar;
 	int _i; _p = _pp; _ppvar = _ppd;
	_cvode_ieq = _ieq;
	for (_i=0; _i < 1; ++_i) {
		_pv[_i] = _pp + _slist1[_i];  _pvdot[_i] = _pp + _dlist1[_i];
		_cvode_abstol(_atollist, _atol, _i);
	}
 }
 
static void _ode_matsol_instance1(_threadargsproto_) {
 _ode_matsol1 (_p, _ppvar, _thread, _nt);
 }
 
static void _ode_matsol(_NrnThread* _nt, _Memb_list* _ml, int _type) {
   double* _p; Datum* _ppvar; Datum* _thread;
   Node* _nd; double _v; int _iml, _cntml;
  _cntml = _ml->_nodecount;
  _thread = _ml->_thread;
  for (_iml = 0; _iml < _cntml; ++_iml) {
    _p = _ml->_data[_iml]; _ppvar = _ml->_pdata[_iml];
    _nd = _ml->_nodelist[_iml];
    v = NODEV(_nd);
 _ode_matsol_instance1(_threadargs_);
 }}

static void initmodel(double* _p, Datum* _ppvar, Datum* _thread, _NrnThread* _nt) {
  int _i; double _save;{
  scount = scount0;
 {
   tspike = - 100000.0 ;
   net_send ( _tqitem, (double*)0, _ppvar[1]._pvoid, t +  0.0 , 1.0 ) ;
   flagOLD = 1.0 ;
   scount = scount0 ;
   d = d0 ;
   p = p0 ;
   alpha_scount = alpha * scount ;
   boltzfactor = exp ( - 1.0 / scounttau ) ;
   if ( boltzman  == 0.0 ) {
     pf = 1.0 / ( 1.0 - boltzfactor ) ;
     }
   else {
     pf = 1.0 ;
     }
   output = 0.0 ;
   }
 
}
}

static void nrn_init(_NrnThread* _nt, _Memb_list* _ml, int _type){
double* _p; Datum* _ppvar; Datum* _thread;
Node *_nd; double _v; int* _ni; int _iml, _cntml;
#if CACHEVEC
    _ni = _ml->_nodeindices;
#endif
_cntml = _ml->_nodecount;
_thread = _ml->_thread;
for (_iml = 0; _iml < _cntml; ++_iml) {
 _p = _ml->_data[_iml]; _ppvar = _ml->_pdata[_iml];
 _tsav = -1e20;
#if CACHEVEC
  if (use_cachevec) {
    _v = VEC_V(_ni[_iml]);
  }else
#endif
  {
    _nd = _ml->_nodelist[_iml];
    _v = NODEV(_nd);
  }
 v = _v;
 initmodel(_p, _ppvar, _thread, _nt);
}
}

static double _nrn_current(double* _p, Datum* _ppvar, Datum* _thread, _NrnThread* _nt, double _v){double _current=0.;v=_v;{
} return _current;
}

static void nrn_cur(_NrnThread* _nt, _Memb_list* _ml, int _type) {
double* _p; Datum* _ppvar; Datum* _thread;
Node *_nd; int* _ni; double _rhs, _v; int _iml, _cntml;
#if CACHEVEC
    _ni = _ml->_nodeindices;
#endif
_cntml = _ml->_nodecount;
_thread = _ml->_thread;
for (_iml = 0; _iml < _cntml; ++_iml) {
 _p = _ml->_data[_iml]; _ppvar = _ml->_pdata[_iml];
#if CACHEVEC
  if (use_cachevec) {
    _v = VEC_V(_ni[_iml]);
  }else
#endif
  {
    _nd = _ml->_nodelist[_iml];
    _v = NODEV(_nd);
  }
 
}
 
}

static void nrn_jacob(_NrnThread* _nt, _Memb_list* _ml, int _type) {
double* _p; Datum* _ppvar; Datum* _thread;
Node *_nd; int* _ni; int _iml, _cntml;
#if CACHEVEC
    _ni = _ml->_nodeindices;
#endif
_cntml = _ml->_nodecount;
_thread = _ml->_thread;
for (_iml = 0; _iml < _cntml; ++_iml) {
 _p = _ml->_data[_iml];
#if CACHEVEC
  if (use_cachevec) {
	VEC_D(_ni[_iml]) += _g;
  }else
#endif
  {
     _nd = _ml->_nodelist[_iml];
	NODED(_nd) += _g;
  }
 
}
 
}

static void nrn_state(_NrnThread* _nt, _Memb_list* _ml, int _type) {
double* _p; Datum* _ppvar; Datum* _thread;
Node *_nd; double _v = 0.0; int* _ni; int _iml, _cntml;
#if CACHEVEC
    _ni = _ml->_nodeindices;
#endif
_cntml = _ml->_nodecount;
_thread = _ml->_thread;
for (_iml = 0; _iml < _cntml; ++_iml) {
 _p = _ml->_data[_iml]; _ppvar = _ml->_pdata[_iml];
 _nd = _ml->_nodelist[_iml];
#if CACHEVEC
  if (use_cachevec) {
    _v = VEC_V(_ni[_iml]);
  }else
#endif
  {
    _nd = _ml->_nodelist[_iml];
    _v = NODEV(_nd);
  }
 v=_v;
{
 {   state(_p, _ppvar, _thread, _nt);
  } {
   scount = scount + output / pf ;
   pf = pf * boltzfactor + 1.0 ;
   alpha_scount = alpha * scount ;
   if ( alpha_scount > 0.0 ) {
     p = p0 / ( bcm_scale * ( alpha_scount - 1.0 ) + 1.0 ) ;
     }
   else {
     p = p ;
     }
   d = d0 * ( bcm_scale * ( alpha_scount - 1.0 ) + 1.0 ) ;
   output = 0.0 ;
   }
}}

}

static void terminal(){}

static void _initlists(){
 double _x; double* _p = &_x;
 int _i; static int _first = 1;
  if (!_first) return;
 _slist1[0] = &(scount) - _p;  _dlist1[0] = &(Dscount) - _p;
_first = 0;
}

#if defined(__cplusplus)
} /* extern "C" */
#endif

#if NMODL_TEXT
static const char* nmodl_filename = "/Users/shkim/Research/neurcomp/Neural-information-processing-and-computations-of-two-input-synapses/Fig2BandFig3/mod_files/BCMthreshold2.mod";
static const char* nmodl_file_text = 
  "COMMENT\n"
  "\n"
  "This is a mechanism for computing sliding BCM threshold based on recent spike count (Benuskova et al. PNAS 2001, Benuskova and Abraham JCNS 2007):\n"
  "\n"
  "alpha_scount = alpha * scount \n"
  "\n"
  "The averaged postsynaptic activity scount expresses the weighted\n"
  "average of the postsynaptic spike count, with the most recent\n"
  "spikes entering the sum with bigger weight than previous\n"
  "ones. \n"
  "\n"
  "output = 1, if there is a postsynaptic spike\n"
  "at a given time, output = 0, otherwise.\n"
  "\n"
  "alpha is the scaling constant.\n"
  "\n"
  "scounttau is the averaging time constant for calculation of alpha_scount.\n"
  "\n"
  "The mechanism should be inserted into soma to calculate the value of alpha_scount and thereby of d and p for all synaptic point processes which use d, p as POINTER variables.\n"
  "At the hoc level d and p have to be set up as POINTER variables to allow the synaptic point process to know the d and p value. Here is an example\n"
  "for setting up POINTER variables for a synaptic object syn: \n"
  "\n"
  "setpointer syn.d, d_BCMthreshold\n"
  "setpointer syn.p, p_BCMthreshold\n"
  "\n"
  "ENDCOMMENT\n"
  "\n"
  "NEURON {\n"
  "	POINT_PROCESS BCMthreshold2\n"
  "	RANGE d0, p0, scount0, scounttau, alpha, alpha_scount, bcm_scale 		:these were GLOBAL before\n"
  "	RANGE d, p, tspike		:these were GLOBAL before\n"
  "}\n"
  "\n"
  "UNITS {\n"
  "	(mV) = (millivolt)\n"
  "}\n"
  "\n"
  "PARAMETER {\n"
  "	\n"
  "	d0		: initial value for the depression factor (additive, non-saturating)\n"
  "	p0		: initial value for the potentiation factor (additive, non-saturating)\n"
  "\n"
  "	scounttau  		: averaging time constant for postsynaptic spike count, e.g. 12000 ms\n"
  "	alpha			: scaling constant\n"
  "	\n"
  "	scount0 		: initial scount = 0 \n"
  "	bcm_scale = 0.2		: scale strength of BCM\n"
  "	boltzman = 0	: initial boltzman = 0 (from Luba)\n"
  "\n"
  "}\n"
  "\n"
  "ASSIGNED {\n"
  "	v (mV)\n"
  "	flagOLD\n"
  "	\n"
  "	alpha_scount	: scount scaled by alpha - sliding (scount-dependent) modification threshold\n"
  "	d			: depression factor (multiplicative to prevent < 0)\n"
  "	p			: potentiation factor (multiplicative)\n"
  "	boltzfactor 	: from Luba\n"
  "	pf				: from Luba\n"
  "	output			: from Luba\n"
  "	tspike (ms)\n"
  "}\n"
  "STATE {\n"
  "	scount			: counter for postsynaptic spikes\n"
  "}\n"
  "INITIAL {\n"
  "	\n"
  "	tspike = -100000\n"
  "	\n"
  ":flag is an implicit argument to NET_RECEIVE. It is an integer, zero by default but a nonzero value can be specified via the second argument to\n"
  ":net_send(). net_send() is used to launch self-events. A self-event comes back to the mechanism that launched it, no NetCon required. \n"
  "\n"
  "	net_send(0, 1)\n"
  "	flagOLD = 1\n"
  "	scount = scount0\n"
  "	d = d0		\n"
  "	p = p0\n"
  "	alpha_scount = alpha*scount\n"
  "	boltzfactor = exp( - 1.0 / scounttau)\n"
  "	if(boltzman == 0) {pf = 1.0 / (1 - boltzfactor)} else {pf = 1.0}\n"
  "	output = 0\n"
  "}\n"
  "\n"
  "BREAKPOINT { \n"
  "	:if (output > 0) {printf(\"entry time=%g output=%g\\n\", t, output)}\n"
  "	scount = scount  + output/pf\n"
  "	SOLVE state METHOD cnexp								:credit to Steffen Platschek\n"
  "  	pf = pf * boltzfactor + 1.0\n"
  "  	alpha_scount = alpha * scount							:scount scaled by alpha\n"
  "	if (alpha_scount > 0) {p = p0/( bcm_scale * (alpha_scount-1) + 1) } else {p = p}\n"
  "	d = d0* (bcm_scale*(alpha_scount-1) +1)\n"
  "    output = 0\n"
  "}\n"
  "\n"
  "DERIVATIVE state{\n"
  "	scount' = -scount / pf\n"
  "}\n"
  "\n"
  ":The items in the argument list in the NET_RECEIVE statement of a synaptic mechanism are actually the elements of a\n"
  ":weight vector. The first element, which is called nc.weight[0],\n"
  ":corresponds to the first item in the NET_RECEIVE argument list--w--which remains constant during a simulation.\n"
  ":The second element is called nc.weight[1], and it corresponds to wE.\n"
  ":The value of this second element DOES change as a consequence of STDP.\n"
  "\n"
  "NET_RECEIVE(w) {		:w is actually not needed for computations\n"
  "	INITIAL {w=0}\n"
  "\n"
  ":When a presynaptic spike occurs, the mechanism receives an event with flag == 0. \n"
  ":The printf statements are purely for diagnostic purposes and can be commented out.\n"
  "\n"
  "	if (flag == 0) {\n"
  "		output = 0 			:no postsynaptic spike\n"
  "	} else if (flag == 2) { 	: postsynaptic spike\n"
  "		tspike = t				: just in case one needs time of the spike\n"
  "		output = 1				:postsynaptic spike \n"
  "		:printf(\"output=%g at time t=tspike=%g\\n\", scount, tspike)\n"
  "	} else { : flag == 1 from INITIAL block\n"
  "\n"
  "		:printf(\"entry flag=%g t=%g\\n\", flag, t)\n"
  "\n"
  ":WATCH (var > thresh) flagvalue is used in a NET_RECEIVE block to specify a condition in the postsynaptic cell\n"
  ":that will generate a self-event with latency 0 and a specified flag value. Generally, WATCH is used to make NEURON\n"
  ":monitor a variable for threshold crossing, and generates a self event with the specified flag value when the threshold\n"
  ":is crossed. If the postsynaptic cell is a biophysical model cell, var is usually local membrane potential (or cai or some\n"
  ":other concentration); if the postsynaptic cell is an artificial spiking cell, var is one of that cell's state variables.\n"
  ":But WATCH could in principle be anything, such as the total number of spikes that a cell has fired, or perhaps even t (time)\n"
  "\n"
  "		WATCH (v > 0) 2 :This mechanism watches postsynaptic membrane potential at the location of the mechanism\n"
  "						:When a postsynaptic spike occurs, the mechanism receives an event with flag == 2\n"
  "	}\n"
  "}\n"
  ;
#endif
