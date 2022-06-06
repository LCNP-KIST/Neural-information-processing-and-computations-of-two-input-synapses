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
 
#define nrn_init _nrn_init__Exp2SynSTDP
#define _nrn_initial _nrn_initial__Exp2SynSTDP
#define nrn_cur _nrn_cur__Exp2SynSTDP
#define _nrn_current _nrn_current__Exp2SynSTDP
#define nrn_jacob _nrn_jacob__Exp2SynSTDP
#define nrn_state _nrn_state__Exp2SynSTDP
#define _net_receive _net_receive__Exp2SynSTDP 
#define state state__Exp2SynSTDP 
 
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
#define tau1 _p[0]
#define tau2 _p[1]
#define e _p[2]
#define dtau _p[3]
#define ptau _p[4]
#define wMax _p[5]
#define start _p[6]
#define wtrack _p[7]
#define d _p[8]
#define p _p[9]
#define i _p[10]
#define g _p[11]
#define tpost _p[12]
#define A _p[13]
#define B _p[14]
#define factor _p[15]
#define presyntime (_p + 16)
#define counter _p[10016]
#define total _p[10017]
#define flagOLD _p[10018]
#define DA _p[10019]
#define DB _p[10020]
#define v _p[10021]
#define _g _p[10022]
#define _tsav _p[10023]
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
 /* some parameters have upper and lower limits */
 static HocParmLimits _hoc_parm_limits[] = {
 "tau2", 1e-09, 1e+09,
 "tau1", 1e-09, 1e+09,
 0,0,0
};
 static HocParmUnits _hoc_parm_units[] = {
 "tau1", "ms",
 "tau2", "ms",
 "e", "mV",
 "dtau", "ms",
 "ptau", "ms",
 "wMax", "uS",
 "start", "ms",
 "A", "uS",
 "B", "uS",
 "i", "nA",
 "g", "uS",
 "tpost", "ms",
 0,0
};
 static double A0 = 0;
 static double B0 = 0;
 static double delta_t = 0.01;
 /* connect global user variables to hoc */
 static DoubScal hoc_scdoub[] = {
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
 
#define _fnc_index 5
 static void _hoc_destroy_pnt(_vptr) void* _vptr; {
   Prop* _prop = ((Point_process*)_vptr)->_prop;
   if (_prop) { _nrn_free_watch(_prop->dparam, 3, 2);}
   if (_prop) { _nrn_free_fornetcon(&(_prop->dparam[_fnc_index]._pvoid));}
   destroy_point_process(_vptr);
}
 
static int _ode_count(int);
static void _ode_map(int, double**, double**, double*, Datum*, double*, int);
static void _ode_spec(_NrnThread*, _Memb_list*, int);
static void _ode_matsol(_NrnThread*, _Memb_list*, int);
 
#define _cvode_ieq _ppvar[6]._i
 static void _ode_matsol_instance1(_threadargsproto_);
 /* connect range variables in _p that hoc is supposed to know about */
 static const char *_mechanism[] = {
 "7.7.0",
"Exp2SynSTDP",
 "tau1",
 "tau2",
 "e",
 "dtau",
 "ptau",
 "wMax",
 "start",
 "wtrack",
 "d",
 "p",
 0,
 "i",
 "g",
 "tpost",
 0,
 "A",
 "B",
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
 	_p = nrn_prop_data_alloc(_mechtype, 10024, _prop);
 	/*initialize range parameters*/
 	tau1 = 0.1;
 	tau2 = 10;
 	e = 0;
 	dtau = 36;
 	ptau = 26;
 	wMax = 0.01;
 	start = 30000;
 	wtrack = 0;
 	d = 0.001;
 	p = 0.003;
  }
 	_prop->param = _p;
 	_prop->param_size = 10024;
  if (!nrn_point_prop_) {
 	_ppvar = nrn_prop_datum_alloc(_mechtype, 7, _prop);
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
 extern int _nrn_netcon_args(void*, double***);
 static void _net_init(Point_process*, double*, double);
 extern Symbol* hoc_lookup(const char*);
extern void _nrn_thread_reg(int, int, void(*)(Datum*));
extern void _nrn_thread_table_reg(int, void(*)(double*, Datum*, Datum*, _NrnThread*, int));
extern void hoc_register_tolerance(int, HocStateTolerance*, Symbol***);
extern void _cvode_abstol( Symbol**, double*, int);

 void _Exp2SynSTDP_reg() {
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
  hoc_register_prop_size(_mechtype, 10024, 7);
  hoc_register_dparam_semantics(_mechtype, 0, "area");
  hoc_register_dparam_semantics(_mechtype, 1, "pntproc");
  hoc_register_dparam_semantics(_mechtype, 2, "netsend");
  hoc_register_dparam_semantics(_mechtype, 3, "watch");
  hoc_register_dparam_semantics(_mechtype, 4, "watch");
  hoc_register_dparam_semantics(_mechtype, 5, "fornetcon");
  hoc_register_dparam_semantics(_mechtype, 6, "cvodeieq");
 	hoc_register_cvode(_mechtype, _ode_count, _ode_map, _ode_spec, _ode_matsol);
 	hoc_register_tolerance(_mechtype, _hoc_state_tol, &_atollist);
 pnt_receive[_mechtype] = _net_receive;
 pnt_receive_init[_mechtype] = _net_init;
 pnt_receive_size[_mechtype] = 4;
 add_nrn_fornetcons(_mechtype, _fnc_index);
 	hoc_register_var(hoc_scdoub, hoc_vdoub, hoc_intfunc);
 	ivoc_help("help ?1 Exp2SynSTDP /Users/shkim/Research/neurcomp/Neural-information-processing-and-computations-of-two-input-synapses/Fig2BandFig3/x86_64/Exp2SynSTDP.mod\n");
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
 static int _slist1[2], _dlist1[2];
 static int state(_threadargsproto_);
 
/*CVODE*/
 static int _ode_spec1 (double* _p, Datum* _ppvar, Datum* _thread, _NrnThread* _nt) {int _reset = 0; {
   DA = - A / tau1 ;
   DB = - B / tau2 ;
   }
 return _reset;
}
 static int _ode_matsol1 (double* _p, Datum* _ppvar, Datum* _thread, _NrnThread* _nt) {
 DA = DA  / (1. - dt*( ( - 1.0 ) / tau1 )) ;
 DB = DB  / (1. - dt*( ( - 1.0 ) / tau2 )) ;
  return 0;
}
 /*END CVODE*/
 static int state (double* _p, Datum* _ppvar, Datum* _thread, _NrnThread* _nt) { {
    A = A + (1. - exp(dt*(( - 1.0 ) / tau1)))*(- ( 0.0 ) / ( ( - 1.0 ) / tau1 ) - A) ;
    B = B + (1. - exp(dt*(( - 1.0 ) / tau2)))*(- ( 0.0 ) / ( ( - 1.0 ) / tau2 ) - B) ;
   }
  return 0;
}
 
static double _watch1_cond(_pnt) Point_process* _pnt; {
 	double* _p; Datum* _ppvar; Datum* _thread; _NrnThread* _nt;
	_thread= (Datum*)0; _nt = (_NrnThread*)_pnt->_vnt;
 	_p = _pnt->_prop->param; _ppvar = _pnt->_prop->dparam;
	v = NODEV(_pnt->node);
	return  ( v ) - ( - 37.0 ) ;
}
 
static void _net_receive (_pnt, _args, _lflag) Point_process* _pnt; double* _args; double _lflag; 
{  double* _p; Datum* _ppvar; Datum* _thread; _NrnThread* _nt;
   int _watch_rm = 0;
   _thread = (Datum*)0; _nt = (_NrnThread*)_pnt->_vnt;   _p = _pnt->_prop->param; _ppvar = _pnt->_prop->dparam;
  if (_tsav > t){ extern char* hoc_object_name(); hoc_execerror(hoc_object_name(_pnt->ob), ":Event arrived out of order. Must call ParallelContext.set_maxstep AFTER assigning minimum NetCon.delay");}
 _tsav = t;   if (_lflag == 1. ) {*(_tqitem) = 0;}
 {
   if ( _lflag  == 0.0 ) {
       if (nrn_netrec_state_adjust && !cvode_active_){
    /* discon state adjustment for cnexp case (rate uses no local variable) */
    double __state = A;
    double __primary = (A + _args[1] * factor) - __state;
     __primary += ( 1. - exp( 0.5*dt*( ( - 1.0 ) / tau1 ) ) )*( - ( 0.0 ) / ( ( - 1.0 ) / tau1 ) - __primary );
    A += __primary;
  } else {
 A = A + _args[1] * factor ;
       }
   if (nrn_netrec_state_adjust && !cvode_active_){
    /* discon state adjustment for cnexp case (rate uses no local variable) */
    double __state = B;
    double __primary = (B + _args[1] * factor) - __state;
     __primary += ( 1. - exp( 0.5*dt*( ( - 1.0 ) / tau2 ) ) )*( - ( 0.0 ) / ( ( - 1.0 ) / tau2 ) - __primary );
    B += __primary;
  } else {
 B = B + _args[1] * factor ;
       }
 _args[2] = t ;
     counter = counter + 1.0 ;
     presyntime [ ((int) counter ) - 1 ] = _args[2] ;
     if ( t > start ) {
       _args[3] = d * exp ( ( tpost - t ) / dtau ) ;
       _args[1] = _args[1] * ( 1.0 - _args[3] ) ;
       if ( _args[1] > 0.0 ) {
         }
       else {
         _args[1] = 0.0 ;
         }
       wtrack = _args[1] ;
       }
     flagOLD = _lflag ;
     }
   else if ( _lflag  == 2.0 ) {
     {int _ifn1, _nfn1; double* _fnargs1, **_fnargslist1;
	_nfn1 = _nrn_netcon_args(_ppvar[_fnc_index]._pvoid, &_fnargslist1);
	for (_ifn1 = 0; _ifn1 < _nfn1; ++_ifn1) {
 	 _fnargs1 = _fnargslist1[_ifn1];
 {
       if ( flagOLD  == _lflag ) {
         }
       else {
         if ( t > start ) {
           {int  _li ;for ( _li = 0 ; _li <= ((int) counter ) - 1 ; _li ++ ) {
             _fnargs1[3] = p * exp ( ( presyntime [ _li ] - t ) / ptau ) ;
             _fnargs1[1] = _fnargs1[1] * ( 1.0 + _fnargs1[3] ) ;
             } }
           if ( _fnargs1[1] < wMax ) {
             }
           else {
             _fnargs1[1] = wMax ;
             wtrack = _fnargs1[1] ;
             }
           {int  _li ;for ( _li = 0 ; _li <= ((int) counter ) - 1 ; _li ++ ) {
             presyntime [ _li ] = - 1e9 ;
             } }
           counter = 0.0 ;
           }
         }
       }
     	}}
 tpost = t ;
     flagOLD = _lflag ;
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
 _args[1] = _args[0] ;
   _args[2] = - 1e9 ;
   _args[3] = 0.0 ;
   wtrack = _args[1] ;
   }
 
static int _ode_count(int _type){ return 2;}
 
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
	for (_i=0; _i < 2; ++_i) {
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
  A = A0;
  B = B0;
 {
   double _ltp ;
 total = 0.0 ;
   if ( tau1 / tau2 > .9999 ) {
     tau1 = .9999 * tau2 ;
     }
   A = 0.0 ;
   B = 0.0 ;
   _ltp = ( tau1 * tau2 ) / ( tau2 - tau1 ) * log ( tau2 / tau1 ) ;
   factor = - exp ( - _ltp / tau1 ) + exp ( - _ltp / tau2 ) ;
   factor = 1.0 / factor ;
   tpost = - 100000.0 ;
   {int  _li ;for ( _li = 0 ; _li <= 10000 - 1 ; _li ++ ) {
     presyntime [ _li ] = - 1e9 ;
     } }
   counter = 0.0 ;
   net_send ( _tqitem, (double*)0, _ppvar[1]._pvoid, t +  0.0 , 1.0 ) ;
   flagOLD = 1.0 ;
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

static double _nrn_current(double* _p, Datum* _ppvar, Datum* _thread, _NrnThread* _nt, double _v){double _current=0.;v=_v;{ {
   g = B - A ;
   i = g * ( v - e ) ;
   }
 _current += i;

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
 _g = _nrn_current(_p, _ppvar, _thread, _nt, _v + .001);
 	{ _rhs = _nrn_current(_p, _ppvar, _thread, _nt, _v);
 	}
 _g = (_g - _rhs)/.001;
 _g *=  1.e2/(_nd_area);
 _rhs *= 1.e2/(_nd_area);
#if CACHEVEC
  if (use_cachevec) {
	VEC_RHS(_ni[_iml]) -= _rhs;
  }else
#endif
  {
	NODERHS(_nd) -= _rhs;
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
  }}}

}

static void terminal(){}

static void _initlists(){
 double _x; double* _p = &_x;
 int _i; static int _first = 1;
  if (!_first) return;
 _slist1[0] = &(A) - _p;  _dlist1[0] = &(DA) - _p;
 _slist1[1] = &(B) - _p;  _dlist1[1] = &(DB) - _p;
_first = 0;
}

#if defined(__cplusplus)
} /* extern "C" */
#endif

#if NMODL_TEXT
static const char* nmodl_filename = "/Users/shkim/Research/neurcomp/Neural-information-processing-and-computations-of-two-input-synapses/Fig2BandFig3/mod_files/Exp2SynSTDP.mod";
static const char* nmodl_file_text = 
  "COMMENT\n"
  "This is a synaptic mechanism implementing a form of stream-specific spike-timing-dependent plasticity (STDP).\n"
  "The mechanism uses simple STDP rule that depends entirely on the interval between the pre- and postsynaptic spikes. This is the nearest-neighbor multiplicative implementation of STDP, that is PRESYNAPTIC CENTRED (Morrison et al. 2008):\n"
  "each presynaptic spike is paired with the last postsynaptic spike and the next postsynaptic spike.\n"
  "\n"
  "When a postsynaptic spike occurs, the mechanism receives an event with flag == 2. The FOR_NETCONS loop iterates over all\n"
  "NetCons that target this particular instance of the synaptic mechanism. It increases each NetCon's weight by\n"
  "a multiplicative factor that depends on the latency between the time of the most recent event (presyn. spike) that was\n"
  "delivered by that NetCon and the time of the postsynaptic spike.\n"
  "\n"
  "When a presynaptic spike occurs, the mechanism receives an event with flag == 0. The weight wE associated with the NetCon\n"
  "that delivered the spike event is depressed by a multiplicative factor that depends on tpost-t, which is the length of time\n"
  "that has elapsed since the most recent postsynaptic spike.\n"
  "\n"
  "The mechanism shown here is a functional model of STDP as opposed to a mechanistic model, which would involve some\n"
  "representation of the biological processes that account for the plasticity. If you need a mechanistic model of STDP that\n"
  "involves processes in the presynaptic terminal, then the equations that describe those processes must be explicitly\n"
  "included in the implementation.\n"
  "\n"
  "This NMODL code also has a WATCH keyword that can be used to monitor membrane potential at the postsynaptic site for\n"
  "occurrence of an action potential (helpful for implementing STDP), and a FOR_NETCONS keyword\n"
  "so that a single instance of an STDP mechanism can deal properly with multiple input streams.\n"
  "\n"
  "FOR_NETCONS(same, args, as, netreceive) { stmt } iterates over all NetCon objects that have this POINT_PROCESS\n"
  "or ARTIFICIAL_CELL as the target. The arglist must have the same number of args as the NET_RECEIVE arg list and must have\n"
  "the same units. If the units are missing, they are inherited by the corresponding arg of the NET_RECEIVE block.\n"
  "If an arg has the same name, it hides the NET_RECEIVE arg. The purpose is to allow straightforward specification of certain\n"
  "kinds of plasticity based on post synaptic state in the context of generalized synapses or artificial cells.\n"
  "\n"
  "ExpSynSTDP is event target, not event source. NetCons convey the fact that an event happened, plus a weight vector and\n"
  "a flag variable. The flag variable is set by code in the target mechanism's NET_RECEIVE block.\n"
  "\n"
  "\n"
  "Two state kinetic scheme synapse described by rise time tau1,\n"
  "and decay time constant tau2. The normalized peak condunductance is 1.\n"
  "Decay time MUST be greater than rise time.\n"
  "\n"
  "The solution of A->G->bath with rate constants 1/tau1 and 1/tau2 is\n"
  " A = a*exp(-t/tau1) and\n"
  " G = a*tau2/(tau2-tau1)*(-exp(-t/tau1) + exp(-t/tau2))\n"
  "	where tau1 < tau2\n"
  "\n"
  "If tau2-tau1 -> 0 then we have a alphasynapse.\n"
  "and if tau1 -> 0 then we have just single exponential decay.\n"
  "\n"
  "The \"factor\" is evaluated in the\n"
  "initial block such that an event of weight 1 generates a\n"
  "peak conductance of 1.\n"
  "\n"
  "Because the solution is a sum of exponentials, the\n"
  "coupled equations can be solved as a pair of independent equations\n"
  "by the more efficient cnexp method.\n"
  "\n"
  "ENDCOMMENT\n"
  "\n"
  "NEURON {\n"
  "	POINT_PROCESS Exp2SynSTDP\n"
  "	RANGE tau1, tau2, e, i, dtau, ptau, wtrack\n"
  "	RANGE wMax 	:hard bound \n"
  "	RANGE d, p	: depression/potentiation factor (additive, non-saturating)\n"
  "	NONSPECIFIC_CURRENT i\n"
  "	RANGE g, tpost, start\n"
  "	:GLOBAL d, p\n"
  "}\n"
  "\n"
  "DEFINE EpspTimes 10000	:to store presynaptic epsp times into an array\n"
  "\n"
  "UNITS {\n"
  "	(nA) = (nanoamp)\n"
  "	(mV) = (millivolt)\n"
  "	(uS) = (microsiemens)\n"
  "}\n"
  "\n"
  "PARAMETER {\n"
  "	tau1=.1 (ms) <1e-9,1e9>\n"
  "	tau2 = 10 (ms) <1e-9,1e9>\n"
  "	e=0	(mV)\n"
  "	dtau = 36 (ms) 			: depression effectiveness time constant (time it takes to fall from 1->0.37)\n"
  "	ptau = 26 (ms) 			: Lin et al. (Eur J Neurosci 2006) dtau 36 ms, ptau 26 ms, Bi & Poo (1998, 2001) dtau 34 ms, ptau 17 ms\n"
  "	wMax = 0.01 (uS)		: hard upper bound \n"
  "	start = 30000 (ms)		: synapses are allowed to start changing after start ms\n"
  "	wtrack = 0\n"
  "	d = 0.001							: depression factor (multiplicative to prevent < 0) \n"
  "	p = 0.003							: potentiation factor (multiplicative) \n"
  "}\n"
  "\n"
  "ASSIGNED {\n"
  "	v (mV)\n"
  "	i (nA)\n"
  "	g (uS)\n"
  "	factor\n"
  "	tpost (ms)\n"
  "	presyntime[EpspTimes] (ms)	: array to store epsp times\n"
  "	counter						: to count epsps until first post-AP comes\n"
  "	total (uS)\n"
  "	flagOLD\n"
  "}\n"
  "\n"
  "STATE {\n"
  "	A (uS)\n"
  "	B (uS)\n"
  "}\n"
  "\n"
  "INITIAL {\n"
  "	LOCAL tp\n"
  "	total = 0\n"
  "	if (tau1/tau2 > .9999) {\n"
  "		tau1 = .9999*tau2\n"
  "	}\n"
  "	A = 0\n"
  "	B = 0\n"
  "	tp = (tau1*tau2)/(tau2 - tau1) * log(tau2/tau1)\n"
  "	factor = -exp(-tp/tau1) + exp(-tp/tau2)\n"
  "	factor = 1/factor\n"
  "	tpost = -100000\n"
  "\n"
  "	FROM i=0 TO EpspTimes-1 {\n"
  "		presyntime[i]=-1e9\n"
  "	}\n"
  "	counter = 0\n"
  "	\n"
  ":flag is an implicit argument to NET_RECEIVE. It is an integer, zero by default but a nonzero value can be specified via the second argument to :net_send(). net_send() is used to launch self-events. A self-event comes back to the mechanism that launched it, no NetCon required. \n"
  "	net_send(0, 1)\n"
  "	flagOLD = 1\n"
  "}\n"
  "\n"
  "BREAKPOINT { \n"
  "	SOLVE state METHOD cnexp\n"
  "	g = B - A\n"
  "	i = g*(v - e)\n"
  "}\n"
  "\n"
  "DERIVATIVE state {\n"
  "	A' = -A/tau1\n"
  "	B' = -B/tau2\n"
  "}\n"
  "\n"
  ":The items in the argument list in the NET_RECEIVE statement of a synaptic mechanism are actually the elements of a\n"
  ":weight vector.  ExpSyn has only one argument in its NET_RECEIVE statement, so a NetCon nc that targets an ExpSyn\n"
  ":has a weight vector with length 1, which may be called nc.weight or nc.weight[0]. A NetCon called nc that targets\n"
  ":an ExpSynSTDP or an Exp2SynSTDP has a weight vector with length 2.  The first element, which is called nc.weight[0],\n"
  ":corresponds to the first item in the NET_RECEIVE argument list--w--which remains constant during a simulation.\n"
  ":The second element is called nc.weight[1], and it corresponds to wE.\n"
  ":The value of this second element DOES change as a consequence of STDP.\n"
  "\n"
  "NET_RECEIVE(w (uS), wE (uS), tpre (ms), X) {\n"
  "	INITIAL { wE = w  tpre = -1e9	 X=0  wtrack=wE}\n"
  "\n"
  ":When a presynaptic spike occurs, the mechanism receives an event with flag == 0. \n"
  ":The printf statements are purely for diagnostic purposes and can be commented out.\n"
  "\n"
  "	if (flag == 0) { : presynaptic spike  (after last post so depress)\n"
  "		:printf(\"entry flag=%g t=%g wInitial=%g wE=%g X=%g tpre=%g tpost=%g\\n\", flag, t, w, wE, X, tpre, tpost)\n"
  "		\n"
  "		A = A + wE*factor\n"
  "		B = B + wE*factor\n"
  "		tpre = t\n"
  "		counter=counter+1				:counting consecutive epsps\n"
  "		presyntime[counter-1]=tpre		:storing epsp time into an array\n"
  "\n"
  ":X*w is the amount by which weight is perturbed from its previous level\n"
  ":\"pre before post\" spiking potentiates, and \"post before pre\" spiking depresses.\n"
  ":Effective synaptic weight is wE \n"
  ":The weight wE associated with the NetCon that delivered the spike \n"
  ":event is depressed by a multiplicative factor that depends on tpost-t, which is \n"
  ":the length of time that has elapsed since the most recent postsynaptic spike.\n"
  "		 \n"
  "		if (t>start) {\n"
  "			X = d*exp((tpost - t)/dtau)\n"
  "			:wE = wE + X*w	:w is initial weight (uS)\n"
  "			wE = wE*(1-X)\n"
  "			if (wE>0) {} else {wE=0}\n"
  "			wtrack = wE	\n"
  "		}\n"
  "		flagOLD = flag\n"
  "	}else if (flag == 2) { : postsynaptic spike\n"
  "		:printf(\"entry flag=%g t=%g tpost=%g\\n\", flag, t, tpost)\n"
  "\n"
  ":The FOR_NETCONS loop iterates over all NetCons that target this particular instance of the synaptic mechanism. It changes each NetCon's X\n"
  ": so that it becomes a potentiation factor depending on the latency between the time of the most recent event (spike) that was delivered by that \n"
  ": NetCon and the time of the postsynaptic spike.\n"
  "\n"
  "		FOR_NETCONS(w1, wE1, tpres, X1) { : also can hide NET_RECEIVE args\n"
  "		:printf(\"entry FOR_NETCONS w1=%g wE1=%g tpres=%g X1=%g\\n\", w1, wE1, tpres, X1)\n"
  "		 \n"
  "		if (flagOLD==flag) {} else {	:for each postsynaptic spike, only 1 presynaptic spike is considered\n"
  "							if (t>start) {\n"
  "								FROM i=0 TO counter-1 {\n"
  "									X1 = p*exp((presyntime[i] - t)/ptau)\n"
  "									wE1 = wE1*(1+X1)\n"
  "								}\n"
  "								:if (presyntime[counter-1]==tpres) {printf(\"Yes the last epsp time is correct-it is %g ms\", presyntime[counter-1])}	\n"
  "								if (wE1<wMax) {} else {wE1=wMax  wtrack=wE1}	:hard bounds \n"
  "								: postsynaptic spike so clear the epsp times array and set counter to 0\n"
  "								FROM i=0 TO counter-1 {\n"
  "									presyntime[i]=-1e9\n"
  "								}\n"
  "								counter = 0	\n"
  "							}\n"
  "						    } \n"
  "		}\n"
  "		tpost = t \n"
  "		:printf(\"scount=%g at time t=tpost=%g\\n\", scount, tpost)\n"
  "		flagOLD=flag\n"
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
  "		WATCH (v > -37) 2 	:This mechanism watches postsynaptic membrane potential at the location of the synapse\n"
  "							:When a postsynaptic spike occurs, the mechanism receives an event with flag == 2\n"
  "	}\n"
  "}\n"
  ;
#endif
