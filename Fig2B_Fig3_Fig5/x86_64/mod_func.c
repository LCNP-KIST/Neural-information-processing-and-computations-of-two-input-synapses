#include <stdio.h>
#include "hocdec.h"
extern int nrnmpi_myid;
extern int nrn_nobanner_;

extern void _BCMthreshold_reg(void);
extern void _BCMthreshold2_reg(void);
extern void _BCMthreshold3_reg(void);
extern void _Exp2SynSTDP_reg(void);
extern void _Exp2SynSTDP_multNNb_globBCM_intscount_precentred_reg(void);
extern void _cad2_reg(void);
extern void _dCaAP_reg(void);
extern void _hh2_reg(void);
extern void _hh3_reg(void);
extern void _it2_reg(void);
extern void _kap_reg(void);
extern void _kca_reg(void);
extern void _kdrca1_reg(void);
extern void _na3_reg(void);
extern void _na3dend_reg(void);
extern void _na3shifted_reg(void);
extern void _sca_reg(void);
extern void _stdp_ca_reg(void);
extern void _stdp_m_reg(void);
extern void _strain_reg(void);

void modl_reg(){
  if (!nrn_nobanner_) if (nrnmpi_myid < 1) {
    fprintf(stderr, "Additional mechanisms from files\n");

    fprintf(stderr," mod_files/BCMthreshold.mod");
    fprintf(stderr," mod_files/BCMthreshold2.mod");
    fprintf(stderr," mod_files/BCMthreshold3.mod");
    fprintf(stderr," mod_files/Exp2SynSTDP.mod");
    fprintf(stderr," mod_files/Exp2SynSTDP_multNNb_globBCM_intscount_precentred.mod");
    fprintf(stderr," mod_files/cad2.mod");
    fprintf(stderr," mod_files/dCaAP.mod");
    fprintf(stderr," mod_files/hh2.mod");
    fprintf(stderr," mod_files/hh3.mod");
    fprintf(stderr," mod_files/it2.mod");
    fprintf(stderr," mod_files/kap.mod");
    fprintf(stderr," mod_files/kca.mod");
    fprintf(stderr," mod_files/kdrca1.mod");
    fprintf(stderr," mod_files/na3.mod");
    fprintf(stderr," mod_files/na3dend.mod");
    fprintf(stderr," mod_files/na3shifted.mod");
    fprintf(stderr," mod_files/sca.mod");
    fprintf(stderr," mod_files/stdp_ca.mod");
    fprintf(stderr," mod_files/stdp_m.mod");
    fprintf(stderr," mod_files/strain.mod");
    fprintf(stderr, "\n");
  }
  _BCMthreshold_reg();
  _BCMthreshold2_reg();
  _BCMthreshold3_reg();
  _Exp2SynSTDP_reg();
  _Exp2SynSTDP_multNNb_globBCM_intscount_precentred_reg();
  _cad2_reg();
  _dCaAP_reg();
  _hh2_reg();
  _hh3_reg();
  _it2_reg();
  _kap_reg();
  _kca_reg();
  _kdrca1_reg();
  _na3_reg();
  _na3dend_reg();
  _na3shifted_reg();
  _sca_reg();
  _stdp_ca_reg();
  _stdp_m_reg();
  _strain_reg();
}
