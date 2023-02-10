#pragma once
#include "utils.hpp"
#include "Containers.hpp"

static const bool debug=true;
static constexpr double func_range_min = 10.;  // Define fitting range
static constexpr double func_range_max = 6500.; // Define fitting range
static constexpr double ptmax_multijet = 1300; // Max pt considered for multijet recoil

static const std::map<TString, TString> input_fnames = {
  {"jes", "rootfiles/jecdataRUN.root"},
};

static const std::map<TString, std::map<TString,TString>> reference_obj_map = {
  {"hjesref",       { {"fname", "jes"}, {"type", "jes"},  {"hname","MODE/etaETAMIN-ETAMAX/herr_l2l3res"},}},
  {"herr",          { {"fname", "jes"}, {"type", "err"},  {"hname","MODE/etaETAMIN-ETAMAX/herr"},}},

  {"herr_ref",      { {"fname", "jes"}, {"type", "err"},  {"hname","MODE/etaETAMIN-ETAMAX/herr_ref"},}},
  {"hrun1",         { {"fname", "jes"}, {"type", "err"},  {"hname","MODE/etaETAMIN-ETAMAX/hrun1"},}},
  {"hjes",          { {"fname", "jes"}, {"type", "err"},  {"hname","MODE/etaETAMIN-ETAMAX/hjes"},}},
  {"herr_l2l3res",  { {"fname", "jes"}, {"type", "err"},  {"hname","MODE/etaETAMIN-ETAMAX/herr_l2l3res"},}},
  {"herr_ref",      { {"fname", "jes"}, {"type", "err"},  {"hname","MODE/etaETAMIN-ETAMAX/herr_ref"},}},
  {"herr_noflv",    { {"fname", "jes"}, {"type", "err"},  {"hname","MODE/etaETAMIN-ETAMAX/herr_noflv"},}},
  {"herr_spr",      { {"fname", "jes"}, {"type", "err"},  {"hname","MODE/etaETAMIN-ETAMAX/herr_spr"},}},
  {"herr_pu",       { {"fname", "jes"}, {"type", "err"},  {"hname","MODE/etaETAMIN-ETAMAX/herr_pu"},}},
};

static const std::map<TString, std::map<TString,TString>> input_hnames_map = {
  {"Resp_zjet_mpf",     { {"fname", "jes"}, {"type", "Resp"}, {"hname","MODE/etaETAMIN-ETAMAX/mpfchs1_zjet_a100"},}},
  {"Resp_zjet_db",      { {"fname", "jes"}, {"type", "Resp"}, {"hname","MODE/etaETAMIN-ETAMAX/ptchs_zjet_a100"},}},

  {"Resp_gamjet_mpf",   { {"fname", "jes"}, {"type", "Resp"}, {"hname","MODE/etaETAMIN-ETAMAX/mpfchs1_gamjet_a100"},}},
  {"Resp_gamjet_db",    { {"fname", "jes"}, {"type", "Resp"}, {"hname","MODE/etaETAMIN-ETAMAX/ptchs_gamjet_a100"},}},

  {"Resp_hadw_mpf",     { {"fname", "jes"}, {"type", "Resp"}, {"hname","MODE/etaETAMIN-ETAMAX/mpfchs1_hadw_a30"},}},
  // {"Resp_hadw_db",      { {"fname", "jes"}, {"type", "Resp"}, {"hname","MODE/etaETAMIN-ETAMAX/ptchs_hadw_a30"},}},
  // {"Resp_hadw_db",      { {"fname", "jes"}, {"type", "Resp"}, {"hname","MODE/etaETAMIN-ETAMAX/mpfchs1_hadw_a30"},}},

  {"Resp_multijet_mpf", { {"fname", "jes"}, {"type", "Resp"}, {"hname","MODE/etaETAMIN-ETAMAX/mpfchs1_multijet_a30"},}},
  {"Resp_multijet_db",  { {"fname", "jes"}, {"type", "Resp"}, {"hname","MODE/etaETAMIN-ETAMAX/ptchs_multijet_a30"},}},

  {"Resp_incljet_mpf",  { {"fname", "jes"}, {"type", "Resp"}, {"hname","MODE/etaETAMIN-ETAMAX/mpfchs1_incjet_a100"},}},
  {"Resp_incljet_db",   { {"fname", "jes"}, {"type", "Resp"}, {"hname","MODE/etaETAMIN-ETAMAX/ptchs_incjet_a100"},}},

};

static const std::map<TString, std::map<TString,TString>> pf_composition_map = {

  {"chf", { {"fname", "jes"}, {"type", "chf"},  {"hname","MODE/etaETAMIN-ETAMAX/chf_cmb_ren"},}},
  {"nhf", { {"fname", "jes"}, {"type", "nhf"},  {"hname","MODE/etaETAMIN-ETAMAX/nhf_cmb_ren"},}},
  {"nef", { {"fname", "jes"}, {"type", "nef"},  {"hname","MODE/etaETAMIN-ETAMAX/nef_cmb_ren"},}},

};

static const std::map<TString, std::map<TString,TString>> recoil_hnames_map = {
  {"Resp_multijet_recoil_mpf", { {"fname", "jes"}, {"type", "Resp"}, {"hname","MODE/etaETAMIN-ETAMAX/crecoil_multijet_a30"},}},
  {"Resp_multijet_recoil_db",  { {"fname", "jes"}, {"type", "Resp"}, {"hname","MODE/etaETAMIN-ETAMAX/crecoil_multijet_a30"},}},
};

static const std::map<TString, std::map<TString,TString>> kfsr_hnames_map = {
  {"fsr_zjet_mpf",     { {"appliesTo", "Resp_zjet_mpf"},     {"fname", "jes"}, {"hname","MODE/etaETAMIN-ETAMAX/fsr/hkfsr3_mpfchs1_zjet"},}},
  {"fsr_gamjet_mpf",   { {"appliesTo", "Resp_gamjet_mpf"},   {"fname", "jes"}, {"hname","MODE/etaETAMIN-ETAMAX/fsr/hkfsr3_mpfchs1_gamjet"},}},
  {"fsr_multijet_mpf", { {"appliesTo", "Resp_multijet_mpf"}, {"fname", "jes"}, {"hname","MODE/etaETAMIN-ETAMAX/fsr/hkfsr3_mpfchs1_multijet"},}},

  {"fsr_zjet_db",      { {"appliesTo", "Resp_zjet_db"},      {"fname", "jes"}, {"hname","MODE/etaETAMIN-ETAMAX/fsr/hkfsr3_ptchs_zjet"},}},
  {"fsr_gamjet_db",    { {"appliesTo", "Resp_gamjet_db"},    {"fname", "jes"}, {"hname","MODE/etaETAMIN-ETAMAX/fsr/hkfsr3_ptchs_gamjet"},}},
  {"fsr_multijet_db",  { {"appliesTo", "Resp_multijet_db"},  {"fname", "jes"}, {"hname","MODE/etaETAMIN-ETAMAX/fsr/hkfsr3_ptchs_multijet"},}},
};

static const std::map<TString, std::map<TString,TString>> sources_hnames_map = {
  {"Resp_uncl_zjet_mpf",        { {"appliesTo", "Resp_zjet_mpf"},     {"fname", "jes"}, {"hname","MODE/etaETAMIN-ETAMAX/fsr/hkfsr3_mpfchs1_zjet_mpfu1"},}},
  {"Resp_add_jet_zjet_mpf",     { {"appliesTo", "Resp_zjet_mpf"},     {"fname", "jes"}, {"hname","MODE/etaETAMIN-ETAMAX/fsr/hkfsr3_mpfchs1_zjet_mpfn1"},}},
  {"Resp_scale_zjet_mpf",       { {"appliesTo", "Resp_zjet_mpf"},     {"fname", "jes"}, {"hname","scale0.2"},}},
  {"Resp_hdmscale_zjet_mpf",    { {"appliesTo", "Resp_zjet_mpf"},     {"fname", "jes"}, {"hname","scale0.2"},}},

  // {"Resp_scale_zjet_db",       { {"appliesTo", "Resp_zjet_db"},     {"fname", "jes"}, {"hname","scale0.2"},}},
  // {"Resp_hdmscale_zjet_db",    { {"appliesTo", "Resp_zjet_db"},     {"fname", "jes"}, {"hname","scale0.2"},}},

  {"Resp_uncl_gamjet_mpf",      { {"appliesTo", "Resp_gamjet_mpf"},   {"fname", "jes"}, {"hname","MODE/etaETAMIN-ETAMAX/fsr/hkfsr3_mpfchs1_gamjet_mpfu1"},}},
  {"Resp_add_jet_gamjet_mpf",   { {"appliesTo", "Resp_gamjet_mpf"},   {"fname", "jes"}, {"hname","MODE/etaETAMIN-ETAMAX/fsr/hkfsr3_mpfchs1_gamjet_mpfn1"},}},
  {"Resp_scale_gamjet_mpf",     { {"appliesTo", "Resp_gamjet_mpf"},   {"fname", "jes"}, {"hname","scale0.5"},}},
  {"Resp_hdmscale_gamjet_mpf",  { {"appliesTo", "Resp_gamjet_mpf"},   {"fname", "jes"}, {"hname","scale0.2"},}},

  {"Resp_scale_gamjet_db",     { {"appliesTo", "Resp_gamjet_db"},   {"fname", "jes"}, {"hname","scale0.5"},}},
  {"Resp_hdmscale_gamjet_db",  { {"appliesTo", "Resp_gamjet_db"},   {"fname", "jes"}, {"hname","scale0.2"},}},

  {"Resp_zee_gamjet_eig0_mpf",  { {"appliesTo", "Resp_gamjet_mpf"},   {"fname", "jes"}, {"hname","MODE/etaETAMIN-ETAMAX/sys/zee_gamjet_eig0"},}},
  {"Resp_zee_gamjet_eig1_mpf",  { {"appliesTo", "Resp_gamjet_mpf"},   {"fname", "jes"}, {"hname","MODE/etaETAMIN-ETAMAX/sys/zee_gamjet_eig1"},}},
  {"Resp_zee_gamjet_eig2_mpf",  { {"appliesTo", "Resp_gamjet_mpf"},   {"fname", "jes"}, {"hname","MODE/etaETAMIN-ETAMAX/sys/zee_gamjet_eig2"},}},

  {"Resp_hadw_fitprob_mpf",     { {"appliesTo", "Resp_hadw_mpf"},     {"fname", "jes"}, {"hname","MODE/etaETAMIN-ETAMAX/sys/hadw_ptave_fitprob"},}},
  {"Resp_hadw_fitprob2_mpf",    { {"appliesTo", "Resp_hadw_mpf"},     {"fname", "jes"}, {"hname","MODE/etaETAMIN-ETAMAX/sys/hadw_ptave_fitprob2"},}},

  {"Resp_uncl_multijet_mpf",    { {"appliesTo", "Resp_multijet_mpf"}, {"fname", "jes"}, {"hname","MODE/etaETAMIN-ETAMAX/fsr/hkfsr3_mpfchs1_multijet_mpfu1"},}},
  {"Resp_add_jet_multijet_mpf", { {"appliesTo", "Resp_multijet_mpf"}, {"fname", "jes"}, {"hname","MODE/etaETAMIN-ETAMAX/fsr/hkfsr3_mpfchs1_multijet_mpfn1"},}},

  {"Resp_uncl_zjet_db",         { {"appliesTo", "Resp_zjet_db"},      {"fname", "jes"}, {"hname","MODE/etaETAMIN-ETAMAX/fsr/hkfsr3_ptchs_zjet_mpfu1"},}},
  {"Resp_add_jet_zjet_db",      { {"appliesTo", "Resp_zjet_db"},      {"fname", "jes"}, {"hname","MODE/etaETAMIN-ETAMAX/fsr/hkfsr3_ptchs_zjet_mpfn1"},}},
  {"Resp_uncl_gamjet_db",       { {"appliesTo", "Resp_gamjet_db"},    {"fname", "jes"}, {"hname","MODE/etaETAMIN-ETAMAX/fsr/hkfsr3_ptchs_gamjet_mpfu1"},}},
  {"Resp_add_jet_gamjet_db",    { {"appliesTo", "Resp_gamjet_db"},    {"fname", "jes"}, {"hname","MODE/etaETAMIN-ETAMAX/fsr/hkfsr3_ptchs_gamjet_mpfn1"},}},
  {"Resp_uncl_multijet_db",     { {"appliesTo", "Resp_multijet_db"},  {"fname", "jes"}, {"hname","MODE/etaETAMIN-ETAMAX/fsr/hkfsr3_ptchs_multijet_mpfu1"},}},
  {"Resp_add_jet_multijet_db",  { {"appliesTo", "Resp_multijet_db"},  {"fname", "jes"}, {"hname","MODE/etaETAMIN-ETAMAX/fsr/hkfsr3_ptchs_multijet_mpfn1"},}},

};

static const std::map<TString, std::map<TString,TString>> shapes_map = {
  {"const",    {{"ispositive", "0"}, {"initial", "0"}, {"freeze", "0"}, {"type","const"},   {"appliesTo","Resp"}, {"form", "1"},}},
  {"ftd",      {{"ispositive", "0"}, {"initial", "0"}, {"freeze", "0"}, {"type","ftd"},     {"appliesTo","Resp"}, {"form", "-0.116-0.6417*pow(x/208.,-0.3051)+23.63/x"},}},
  {"fp",       {{"ispositive", "0"}, {"initial", "0"}, {"freeze", "0"}, {"type","fp"},      {"appliesTo","Resp"}, {"form", "-0.8295"},}},
  {"fhx",      {{"ispositive", "0"}, {"initial", "0"}, {"freeze", "0"}, {"type","fhx"},     {"appliesTo","Resp"}, {"form", "0.8904+1.082*pow(x/1408,1.204)/(1+pow(x/1408,1.204))*(1-pow(x/1408,-1.204))"},}},
  {"fhh",      {{"ispositive", "0"}, {"initial", "0"}, {"freeze", "0"}, {"type","fhh"},     {"appliesTo","Resp"}, {"form", "-0.7938-0.5798*pow(x/396.1,1.412)/(1+pow(x/396.1,1.412))*(1-pow(x/396.1,-1.412))"},}},
  {"feh",      {{"ispositive", "0"}, {"initial", "0"}, {"freeze", "0"}, {"type","feh"},     {"appliesTo","Resp"}, {"form", "-0.2603-0.2196*pow(x/409.4,1.276)/(1+pow(x/409.4,1.276))*(1-pow(x/409.4,-1.276))"},}},
  {"fhw",      {{"ispositive", "0"}, {"initial", "0"}, {"freeze", "0"}, {"type","fhw"},     {"appliesTo","Resp"}, {"form", "0.3*(0.9526-0.3883*(1+(pow(x/1285,2.46)-1)/(pow(x/1285,2.46)+1))+18.1/x-2.062*log(x)/x)"},}},
  {"fl1",      {{"ispositive", "0"}, {"initial", "0"}, {"freeze", "0"}, {"type","fl1"},     {"appliesTo","Resp"}, {"form", "100*(1-(0.350077+0.553560*log(x)-0.0527681*pow(log(x),2))/x-1)"},}},
  {"ftd-ftm",  {{"ispositive", "0"}, {"initial", "0"}, {"freeze", "0"}, {"type","ftd-ftm"}, {"appliesTo","Resp"}, {"form", "3*((-0.116-0.6417*pow(x/208.,-0.3051)+23.63/x)-(0.2683-0.6994*pow(x/208.,-0.3051)+18.49/x))"},}},
  {"f1q3-1",   {{"ispositive", "0"}, {"initial", "0"}, {"freeze", "0"}, {"type","f1q3-1"},  {"appliesTo","Resp"}, {"form", "0.01*(0.7966+0.9311*(pow(0.01*x,-1)-1))"},}},
};

static const std::map<TString, std::map<TString,TString>> shapes_pf_map = {
  {"ftd_chf",  {{"ispositive", "0"}, {"initial", "0"}, {"freeze", "0"}, {"type","ftd"}, {"appliesTo","chf"}, {"form", "1.982-2.678*(1+(pow(x/47.02,0.262)-1)/(pow(x/47.02,0.262)+1))+0.1494*pow(x,+0.3)-3.097/x"},}},
  {"ftd_nhf",  {{"ispositive", "0"}, {"initial", "0"}, {"freeze", "0"}, {"type","ftd"}, {"appliesTo","nhf"}, {"form", "-0.01022-0.1962*(1+(pow(x/4000,3.071)-1)/(pow(x/4000,3.071)+1))+0.04211*pow(x,+0.3)+0.01005/x"},}},
  {"ftd_nef",  {{"ispositive", "0"}, {"initial", "0"}, {"freeze", "0"}, {"type","ftd"}, {"appliesTo","nef"}, {"form", "0.07453+0.1457*(1+(pow(x/1131,-3.68)-1)/(pow(x/1131,-3.68)+1))-0.4155*pow(x,-0.3)-1.878/x"},}},

  {"fp_chf",   {{"ispositive", "0"}, {"initial", "0"}, {"freeze", "0"}, {"type","fp"}, {"appliesTo","chf"}, {"form", "0.3333+0.7433*(1+(pow(x/1023,0.3926)-1)/(pow(x/1023,0.3926)+1))-0.09446*pow(x,0.2883)"},}},
  {"fp_nhf",   {{"ispositive", "0"}, {"initial", "0"}, {"freeze", "0"}, {"type","fp"}, {"appliesTo","nhf"}, {"form", "0.07395+0*(1+(pow(x/1000,258.2)-1)/(pow(x/1000,258.2)+1))+1.223e-05*pow(x,1.158)"},}},
  {"fp_nef",   {{"ispositive", "0"}, {"initial", "0"}, {"freeze", "0"}, {"type","fp"}, {"appliesTo","nef"}, {"form", "2.283+0*(1+(pow(x/1000,1.302)-1)/(pow(x/1000,1.302)+1))-2.738*pow(x,0.002452)"},}},

  {"fhx_chf",  {{"ispositive", "0"}, {"initial", "0"}, {"freeze", "0"}, {"type","fhx"}, {"appliesTo","chf"}, {"form", "-0.0637-0.2811*(1+(pow(x/4531,-0.3172)-1)/(pow(x/4531,-0.3172)+1))+1.071*pow(x,-0.153)"},}},
  {"fhx_nhf",  {{"ispositive", "0"}, {"initial", "0"}, {"freeze", "0"}, {"type","fhx"}, {"appliesTo","nhf"}, {"form", "-0.295+0.09444*(1+(pow(x/2713,0.06437)-1)/(pow(x/2713,0.06437)+1))+[p4]*pow(x,0.2845)"},}},
  {"fhx_nef",  {{"ispositive", "0"}, {"initial", "0"}, {"freeze", "0"}, {"type","fhx"}, {"appliesTo","nef"}, {"form", "0.05474-0.003141*(1+(pow(x/798.6,78.84)-1)/(pow(x/798.6,78.84)+1))-0.000957*pow(x,0.76)"},}},

  {"fhh_chf",  {{"ispositive", "0"}, {"initial", "0"}, {"freeze", "0"}, {"type","fhh"}, {"appliesTo","chf"}, {"form", "0.1552-0.04221*(1+(pow(x/315,2.787)-1)/(pow(x/315,2.787)+1))-0.06628*pow(x,-0.2572)"},}},
  {"fhh_nhf",  {{"ispositive", "0"}, {"initial", "0"}, {"freeze", "0"}, {"type","fhh"}, {"appliesTo","nhf"}, {"form", "-0.2746-0.6358*(1+(pow(x/9664,0.6547)-1)/(pow(x/9664,0.6547)+1))+0.05559*pow(x,0.1816)"},}},
  {"fhh_nef",  {{"ispositive", "0"}, {"initial", "0"}, {"freeze", "0"}, {"type","fhh"}, {"appliesTo","nef"}, {"form", "0.4158-2.14*(1+(pow(x/9426,0.1723)-1)/(pow(x/9426,0.1723)+1))+0.4111*pow(x,0.1937)"},}},

  {"feh_chf",  {{"ispositive", "0"}, {"initial", "0"}, {"freeze", "0"}, {"type","feh"}, {"appliesTo","chf"}, {"form", "0.06085+0*(1+(pow(x/1000,1.3)-1)/(pow(x/1000,1.3)+1))-0.008137*pow(x,0.2135)"},}},
  {"feh_nhf",  {{"ispositive", "0"}, {"initial", "0"}, {"freeze", "0"}, {"type","feh"}, {"appliesTo","nhf"}, {"form", "-0.03458+0*(1+(pow(x/1713,274.8)-1)/(pow(x/1713,274.8)+1))+0.01665*pow(x,0.2426)"},}},
  {"feh_nef",  {{"ispositive", "0"}, {"initial", "0"}, {"freeze", "0"}, {"type","feh"}, {"appliesTo","nef"}, {"form", "-0.02364+0*(1+(pow(x/1481,246.2)-1)/(pow(x/1481,246.2)+1))-0.009737*pow(x,0.2576)"},}},

  {"fhw_chf",  {{"ispositive", "0"}, {"initial", "0"}, {"freeze", "0"}, {"type","fhw"}, {"appliesTo","chf"}, {"form", "-0.2176+1.064e-05*pow(x,1.373)+0/x"},}},
  {"fhw_nhf",  {{"ispositive", "0"}, {"initial", "0"}, {"freeze", "0"}, {"type","fhw"}, {"appliesTo","nhf"}, {"form", "-5.151+4.495*pow(x,0.03335)-12.3/x"},}},
  {"fhw_nef",  {{"ispositive", "0"}, {"initial", "0"}, {"freeze", "0"}, {"type","fhw"}, {"appliesTo","nef"}, {"form", "0.8417-0.2605*pow(x,0.2289)+2.426/x"},}},

};
