/*!
 * \file   ZMATNorton.cxx
 * \brief  This file declares the ZMAT interface for the Norton behaviour law
 * \author Helfer Thomas
 * \date   23 / 11 / 06
 */


#include<string>
#include<vector>
#include<algorithm>

#include"External_parameter.h"
#include"Verbose.h"
#include"MFront/ZMAT/ZMATUndefs.hxx"

#include"TFEL/Material/Norton.hxx"
#include"MFront/ZMAT/ZMATNorton.hxx"

Z_START_NAMESPACE;

ZMATNorton::ZMATNorton()
  : obp(tfel::material::None)
{
  this->zbb_keep_ep = &this->local_ep_list;
#ifdef _WIN64
  ZMAT_GLOBAL_STORAGE::zmat_once();
  ZMAT_GLOBAL_STORAGE& zmat_globals = *thread_zmat_globals;
#endif
} // end of ZMATNorton::ZMATNorton()

void ZMATNorton::initialize(ASCII_FILE& file,int dim,LOCAL_INTEGRATION* integ){
  BEHAVIOR::initialize(file, dim,integ);
  using namespace std;
  int keep_verbose = ZSET::Verbose;
  this->coord.resize(dim);
  this->sig.initialize(this,"sig",this->tsz(),1);
  this->eto.initialize(this,"eto",this->tsz(),1);
  this->tg_mat.resize(this->tsz(), this->tsz());
  // initialisation dispatch
  if(this->tsz()==6){
    this->mprops.resize(4);
    this->eel.initialize(this,"ElasticStrain",this->tsz(),1);
    this->p.initialize(this,"EquivalentViscoplasticStrain",1,1);
    for(;;){
      STRING str=file.getSTRING();
      if(this->base_read(str,file)){
      } else if((str=="**model_coef")||(str=="**material_properties")){
	this->initializeMaterialProperties3D(file);
      } else if(str=="**parameters"){
	this->initializeParameters3D(file);
      } else if(str=="**out_of_bounds_policy"){
	STRING policy=file.getSTRING();
	if(policy=="None"){
	  this->obp=tfel::material::None;
	} else if(policy=="Strict"){
	  this->obp=tfel::material::Strict;
	} else if(policy=="Warning"){
	  this->obp=tfel::material::Warning;
	} else {
	  INPUT_ERROR("unknown policy '"+policy+"'");
	}
      } else if (str.start_with("***")){
	file.back();
	break;
      } else {
	INPUT_ERROR("Invalid keyword '"+str+"'");
      }
    }
  } else if(this->tsz()==4){
    this->mprops.resize(4);
    this->eel.initialize(this,"ElasticStrain",this->tsz(),1);
    this->p.initialize(this,"EquivalentViscoplasticStrain",1,1);
    for(;;){
      STRING str=file.getSTRING();
      if(this->base_read(str,file)){
      } else if((str=="**model_coef")||(str=="**material_properties")){
	this->initializeMaterialProperties2D(file);
      } else if(str=="**parameters"){
	this->initializeParameters2D(file);
      } else if(str=="**out_of_bounds_policy"){
	STRING policy=file.getSTRING();
	if(policy=="None"){
	  this->obp=tfel::material::None;
	} else if(policy=="Strict"){
	  this->obp=tfel::material::Strict;
	} else if(policy=="Warning"){
	  this->obp=tfel::material::Warning;
	} else {
	  INPUT_ERROR("unknown policy '"+policy+"'");
	}
      } else if (str.start_with("***")){
	file.back();
	break;
      } else {
	INPUT_ERROR("Invalid keyword '"+str+"'");
      }
    }
  } else if(this->tsz()==3){
    this->mprops.resize(4);
    this->eel.initialize(this,"ElasticStrain",this->tsz(),1);
    this->p.initialize(this,"EquivalentViscoplasticStrain",1,1);
    for(;;){
      STRING str=file.getSTRING();
      if(this->base_read(str,file)){
      } else if((str=="**model_coef")||(str=="**material_properties")){
	this->initializeMaterialProperties1D(file);
      } else if(str=="**parameters"){
	this->initializeParameters1D(file);
      } else if(str=="**out_of_bounds_policy"){
	STRING policy=file.getSTRING();
	if(policy=="None"){
	  this->obp=tfel::material::None;
	} else if(policy=="Strict"){
	  this->obp=tfel::material::Strict;
	} else if(policy=="Warning"){
	  this->obp=tfel::material::Warning;
	} else {
	  INPUT_ERROR("unknown policy '"+policy+"'");
	}
      } else if (str.start_with("***")){
	file.back();
	break;
      } else {
	INPUT_ERROR("Invalid keyword '"+str+"'");
      }
    }
  } else {
    ERROR("Invalid tensor size");
  }
  this->temperature_position = EXTERNAL_PARAM::rank_of_nodal_ip("temperature");
  if(this->temperature_position==-1){;
    INPUT_ERROR("temperature is not defined");
  }
  // check that all material properties were initialised
  for(int pc=0;pc!=this->mprops.size();++pc){
    if(!this->mprops[pc].ok()){
      INPUT_ERROR("Some material properties were not initialised");
    }
  }
  ZSET::Verbose = keep_verbose;
} // end of ZMATNorton::initialize

void
ZMATNorton::initializeMaterialProperties3D(ASCII_FILE& file){
  using std::find;
  const STRING all_mp_names[4] = {
    "NortonCoefficient",
    "NortonExponent",
    "PoissonRatio",
    "YoungModulus"};
  const STRING mp_names3D[4] = {"YoungModulus",
				"PoissonRatio",
				"NortonCoefficient",
				"NortonExponent"};
  for(;;){
    STRING str=file.getSTRING();
    if(str[0]=='*'){
      file.back();
      break;
    }
    if(find(all_mp_names,all_mp_names+4,str)==all_mp_names+4){
      INPUT_ERROR("No material property named '"+str+"'");
    }
    const STRING * const pmat = find(mp_names3D,mp_names3D+4,str);
    if(this->mprops[pmat-mp_names3D].ok()){
      INPUT_ERROR("material property '"+str+"' already defined");
    }
    this->mprops[pmat-mp_names3D].read(str,file,this);
  }
} // end of ZMATNorton::initializeMaterialProperties3D

void
ZMATNorton::initializeMaterialProperties2D(ASCII_FILE& file){
  using std::find;
  const STRING all_mp_names[4] = {
    "NortonCoefficient",
    "NortonExponent",
    "PoissonRatio",
    "YoungModulus"};
  const STRING mp_names2D[4] = {"YoungModulus",
				"PoissonRatio",
				"NortonCoefficient",
				"NortonExponent"};
  for(;;){
    STRING str=file.getSTRING();
    if(str[0]=='*'){
      file.back();
      break;
    }
    if(find(all_mp_names,all_mp_names+4,str)==all_mp_names+4){
      INPUT_ERROR("No material property named '"+str+"'");
    }
    const STRING * const pmat = find(mp_names2D,mp_names2D+4,str);
    if(this->mprops[pmat-mp_names2D].ok()){
      INPUT_ERROR("material property '"+str+"' already defined");
    }
    this->mprops[pmat-mp_names2D].read(str,file,this);
  }
} // end of ZMATNorton::initializeMaterialProperties2D

void
ZMATNorton::initializeMaterialProperties1D(ASCII_FILE& file){
  using std::find;
  const STRING all_mp_names[4] = {
    "NortonCoefficient",
    "NortonExponent",
    "PoissonRatio",
    "YoungModulus"};
  const STRING mp_names1D[4] = {"YoungModulus",
				"PoissonRatio",
				"NortonCoefficient",
				"NortonExponent"};
  for(;;){
    STRING str=file.getSTRING();
    if(str[0]=='*'){
      file.back();
      break;
    }
    if(find(all_mp_names,all_mp_names+4,str)==all_mp_names+4){
      INPUT_ERROR("No material property named '"+str+"'");
    }
    const STRING * const pmat = find(mp_names1D,mp_names1D+4,str);
    if(this->mprops[pmat-mp_names1D].ok()){
      INPUT_ERROR("material property '"+str+"' already defined");
    }
    this->mprops[pmat-mp_names1D].read(str,file,this);
  }
} // end of ZMATNorton::initializeMaterialProperties1D

void
ZMATNorton::initializeParameters3D(ASCII_FILE& file){
  for(;;){
    STRING str=file.getSTRING();
    if(str[0]=='*'){
      file.back();
      break;
    } else if(str=="minimal_time_step_scaling_factor"){
      const double value=file.getdouble();
      tfel::material::NortonParametersInitializer::get().minimal_time_step_scaling_factor = value;
    } else if(str=="maximal_time_step_scaling_factor"){
      const double value=file.getdouble();
      tfel::material::NortonParametersInitializer::get().maximal_time_step_scaling_factor = value;
    } else if(str=="theta"){
      const double value=file.getdouble();
      tfel::material::NortonParametersInitializer::get().theta = value;
    } else if(str=="epsilon"){
      const double value=file.getdouble();
      tfel::material::NortonParametersInitializer::get().epsilon = value;
    } else if(str=="iterMax"){
      const unsigned short value=static_cast<unsigned short>(file.getint());
      tfel::material::NortonParametersInitializer::get().iterMax = value;
    } else {
      INPUT_ERROR("invalid parameter name '"+str+"'");
    }
  }
}

void
ZMATNorton::initializeParameters2D(ASCII_FILE& file){
  for(;;){
    STRING str=file.getSTRING();
    if(str[0]=='*'){
      file.back();
      break;
    } else if(str=="minimal_time_step_scaling_factor"){
      const double value=file.getdouble();
      tfel::material::NortonParametersInitializer::get().minimal_time_step_scaling_factor = value;
    } else if(str=="maximal_time_step_scaling_factor"){
      const double value=file.getdouble();
      tfel::material::NortonParametersInitializer::get().maximal_time_step_scaling_factor = value;
    } else if(str=="theta"){
      const double value=file.getdouble();
      tfel::material::NortonParametersInitializer::get().theta = value;
    } else if(str=="epsilon"){
      const double value=file.getdouble();
      tfel::material::NortonParametersInitializer::get().epsilon = value;
    } else if(str=="iterMax"){
      const unsigned short value=static_cast<unsigned short>(file.getint());
      tfel::material::NortonParametersInitializer::get().iterMax = value;
    } else {
      INPUT_ERROR("invalid parameter name '"+str+"'");
    }
  }
}

void
ZMATNorton::initializeParameters1D(ASCII_FILE& file){
  for(;;){
    STRING str=file.getSTRING();
    if(str[0]=='*'){
      file.back();
      break;
    } else if(str=="minimal_time_step_scaling_factor"){
      const double value=file.getdouble();
      tfel::material::NortonParametersInitializer::get().minimal_time_step_scaling_factor = value;
    } else if(str=="maximal_time_step_scaling_factor"){
      const double value=file.getdouble();
      tfel::material::NortonParametersInitializer::get().maximal_time_step_scaling_factor = value;
    } else if(str=="theta"){
      const double value=file.getdouble();
      tfel::material::NortonParametersInitializer::get().theta = value;
    } else if(str=="epsilon"){
      const double value=file.getdouble();
      tfel::material::NortonParametersInitializer::get().epsilon = value;
    } else if(str=="iterMax"){
      const unsigned short value=static_cast<unsigned short>(file.getint());
      tfel::material::NortonParametersInitializer::get().iterMax = value;
    } else {
      INPUT_ERROR("invalid parameter name '"+str+"'");
    }
  }
}

INTEGRATION_RESULT* ZMATNorton::integrate(MAT_DATA& mdat,const VECTOR& delta_grad,
					  MATRIX*& tg_matrix,int flags){
  int keep_verbose  = ZSET::Verbose; 
  CLOCK* keep_clock = ZSET::stored_thread_zbase_globals->ptr()->active_clock; 
  tg_matrix = &(this->tg_mat);
  this->set_var_aux_to_var_aux_ini();
  this->set_var_int_to_var_int_ini();
  LIST<EXTERNAL_PARAM*>* ep_save = &EXTERNAL_PARAM::Get_EP_list();
  EXTERNAL_PARAM::set_EP_list(this->zbb_keep_ep);
  if(!this->curr_ext_param){
    this->curr_ext_param = *mdat.param_set();
  }
  this->calc_local_coefs();
  INTEGRATION_RESULT* r = NULL;
  try{
    if(this->tsz()==6){
      r=this->callMFrontBehaviour3D(mdat,delta_grad,tg_matrix,flags);
    } else if(this->tsz()==4){
      r=this->callMFrontBehaviour2D(mdat,delta_grad,tg_matrix,flags);
    } else if(this->tsz()==3){
      r=this->callMFrontBehaviour1D(mdat,delta_grad,tg_matrix,flags);
    } else {
      ERROR("Invalid tensor size");
    }
  }
  catch(std::exception& e){
    static INTEGRATION_RESULT bad_result;
    *(*ZSET::stored_thread_zbase_globals).ptr()->out_msg << e.what() << endl;
    bad_result.set_error(INTEGRATION_RESULT::UNDEFINED_BEHAVIOR);
    ZSET::Verbose = keep_verbose; 
    ZSET::stored_thread_zbase_globals->ptr()->active_clock = keep_clock; 
    return &bad_result;
  }
  catch(...){
    static INTEGRATION_RESULT bad_result;
    bad_result.set_error(INTEGRATION_RESULT::UNDEFINED_BEHAVIOR);
    ZSET::Verbose = keep_verbose; 
    ZSET::stored_thread_zbase_globals->ptr()->active_clock = keep_clock; 
    return &bad_result;
  }
  if(r!=NULL){
    ZSET::Verbose = keep_verbose; 
    ZSET::stored_thread_zbase_globals->ptr()->active_clock = keep_clock; 
    return r;
  }
  this->update_var_aux();
  this->zbb_keep_ep = &EXTERNAL_PARAM::Get_EP_list();
  EXTERNAL_PARAM::set_EP_list(ep_save);
  ZSET::Verbose = keep_verbose;
  ZSET::stored_thread_zbase_globals->ptr()->active_clock = keep_clock;
  return NULL;
} // end of ZMATNorton::integrate

INTEGRATION_RESULT* ZMATNorton::callMFrontBehaviour3D(MAT_DATA& mdat,const VECTOR& delta_grad,
						      MATRIX*& tg_matrix,int flags){
  typedef tfel::material::Norton<tfel::material::ModellingHypothesis::TRIDIMENSIONAL,double,false> Norton;

  ....
  auto smflag = Norton::STANDARDTANGENTOPERATOR;
  auto smtype = Norton::NOSTIFFNESSREQUESTED;
  if(flags&CALC_TG_MATRIX){
    smtype = Norton::CONSISTENTTANGENTOPERATOR;
  }
  Norton b(this->sig,stran,dstran,this->mprops,mdat,this->temperature_position,
	   this->evs_positions,ZSET::stored_thread_zbase_globals->ptr()->active_clock->get_dtime());
  b.initialize();
  if(b.integrate(smflag,smtype)!=Norton::SUCCESS){
    static INTEGRATION_RESULT bad_result;
    bad_result.set_error(INTEGRATION_RESULT::UNDEFINED_BEHAVIOR);
    return &bad_result;
  }
  b.ZMATexportStateData(this->sig,mdat);
  if(smtype!=Norton::NOSTIFFNESSREQUESTED){
    zmat::ZMATInterface::convert(*tg_matrix,b.getTangentOperator());
  }
  return nullptr;
} // end of ZMATNorton::callMFrontBehaviour3D

INTEGRATION_RESULT*
ZMATNorton::callMFrontBehaviour2D(MAT_DATA& mdat,
				  const VECTOR& delta_grad,
				  MATRIX*& tg_matrix,
				  int flags){
  typedef tfel::material::ModellingHypothesis ModellingHypothesis;
  typedef tfel::material::Norton<ModellingHypothesis::GENERALISEDPLANESTRAIN,double,false> Norton;
  using tfel::math::st2tost2;
  // strain and strain increment
  double stran[4];
  double dstran[4];
  stran[0] = this->eto[0]-delta_grad[0];
  stran[1] = this->eto[1]-delta_grad[1];
  stran[2] = this->eto[2]-delta_grad[2];
  stran[3] = this->eto[3]-delta_grad[3];
  dstran[0] = delta_grad[0];
  dstran[1] = delta_grad[1];
  dstran[2] = delta_grad[2];
  dstran[3] = delta_grad[3];
  Norton::SMFlag smflag = Norton::STANDARDTANGENTOPERATOR;
  // tangent operator type
  Norton::SMType smtype = Norton::NOSTIFFNESSREQUESTED;
  if(flags&CALC_TG_MATRIX){
    smtype = Norton::CONSISTENTTANGENTOPERATOR;
  }
  Norton b(this->sig,stran,dstran,this->mprops,mdat,this->temperature_position,
	   this->evs_positions,ZSET::stored_thread_zbase_globals->ptr()->active_clock->get_dtime());
  b.initialize();
  if(b.integrate(smflag,smtype)!=Norton::SUCCESS){
    static INTEGRATION_RESULT bad_result;
    bad_result.set_error(INTEGRATION_RESULT::UNDEFINED_BEHAVIOR);
    return &bad_result;
  }
  b.ZMATexportStateData(this->sig,mdat);
  if(smtype!=Norton::NOSTIFFNESSREQUESTED){
    zmat::ZMATInterface::convert(*tg_matrix,b.getTangentOperator());
  }
  return nullptr;
} // end of ZMATNorton::callMFrontBehaviour2D

INTEGRATION_RESULT*
ZMATNorton::callMFrontBehaviour1D(MAT_DATA& mdat,
				  const VECTOR& delta_grad,
				  MATRIX*& tg_matrix,
				  int flags){
  typedef tfel::material::ModellingHypothesis ModellingHypothesis;
  typedef tfel::material::Norton<ModellingHypothesis::AXISYMMETRICALGENERALISEDPLANESTRAIN,double,false> Norton;
  using tfel::math::st2tost2;
  // strain and strain increment
  double stran[3];
  double dstran[3];
  stran[0] = this->eto[0]-delta_grad[0];
  stran[1] = this->eto[1]-delta_grad[1];
  stran[2] = this->eto[2]-delta_grad[2];
  dstran[0] = delta_grad[0];
  dstran[1] = delta_grad[1];
  dstran[2] = delta_grad[2];
  Norton::SMFlag smflag = Norton::STANDARDTANGENTOPERATOR;
  // tangent operator type
  Norton::SMType smtype = Norton::NOSTIFFNESSREQUESTED;
  if(flags&CALC_TG_MATRIX){
    smtype = Norton::CONSISTENTTANGENTOPERATOR;
  }
  Norton b(this->sig,stran,dstran,this->mprops,mdat,this->temperature_position,
	   this->evs_positions,ZSET::stored_thread_zbase_globals->ptr()->active_clock->get_dtime());
  b.initialize();
  if(b.integrate(smflag,smtype)!=Norton::SUCCESS){
    static INTEGRATION_RESULT bad_result;
    bad_result.set_error(INTEGRATION_RESULT::UNDEFINED_BEHAVIOR);
    return &bad_result;
  }
  b.ZMATexportStateData(this->sig,mdat);
  if(smtype!=Norton::NOSTIFFNESSREQUESTED){
    zmat::ZMATInterface::convert(*tg_matrix,b.getTangentOperator());
  }
  return nullptr;
} // end of ZMATNorton::callMFrontBehaviour1D


ZMATNorton::~ZMATNorton(){
} // end of ZMATNorton::~ZMATNorton

BEHAVIOR_READER(ZMATNorton,Norton)

Z_END_NAMESPACE;
