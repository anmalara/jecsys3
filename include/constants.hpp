#pragma once
#include "utils.hpp"
#include "Containers.hpp"

static const bool debug=true;
static constexpr double func_range_min = 10.;  // Define fitting range
static constexpr double func_range_max = 6500.; // Define fitting range
static constexpr double ptmin_multijet = -1; // Min pt considered for multijet recoil
static constexpr double ptmax_multijet = 1300; // Max pt considered for multijet recoil

static constexpr double ptmin_pf = 40; // Min pt considered for pf composition
static constexpr double ptmax_pf = 1000; // Max pt considered for pf composition

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

  {"chf",               { {"fname", "jes"}, {"type", "chf"},  {"hname","MODE/etaETAMIN-ETAMAX/chf_cmb_ren"},}},
  {"nhf",               { {"fname", "jes"}, {"type", "nhf"},  {"hname","MODE/etaETAMIN-ETAMAX/nhf_cmb_ren"},}},
  {"nef",               { {"fname", "jes"}, {"type", "nef"},  {"hname","MODE/etaETAMIN-ETAMAX/nef_cmb_ren"},}},

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

static const std::map<TString, std::map<TString,TString>> nuisances_map = {
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
  {"fl1",      {{"ispositive", "0"}, {"initial", "0"}, {"freeze", "0"}, {"type","fl1"},     {"appliesTo","Resp"}, {"form", "(1-(0.350077+0.553560*log(x)-0.0527681*pow(log(x),2))/x-1)"},}},
  {"ftd-ftm",  {{"ispositive", "0"}, {"initial", "0"}, {"freeze", "0"}, {"type","ftd-ftm"}, {"appliesTo","Resp"}, {"form", "3*((-0.116-0.6417*pow(x/208.,-0.3051)+23.63/x)-(0.2683-0.6994*pow(x/208.,-0.3051)+18.49/x))"},}},
  {"f1q3-1",   {{"ispositive", "0"}, {"initial", "0"}, {"freeze", "0"}, {"type","f1q3-1"},  {"appliesTo","Resp"}, {"form", "0.01*(0.7966+0.9311*(pow(0.01*x,-1)-1))"},}},

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


  {"hadHcalp3",              {{"ispositive", "0"}, {"initial", "0"}, {"freeze", "0"}, {"type","hadHcalp3"},          {"appliesTo","Resp"}, {"form", "+2.792e+00-2.377e+00*log(x)+6.981e-01*pow(log(x),2)-7.477e-02*pow(log(x),3)+2.779e-03*pow(log(x),4)"},}},
  {"hadHcalp3_chf",          {{"ispositive", "0"}, {"initial", "0"}, {"freeze", "0"}, {"type","hadHcalp3"},          {"appliesTo","chf"},  {"form", "-2.368e+00+1.749e+00*log(x)-4.589e-01*pow(log(x),2)+4.307e-02*pow(log(x),3)-1.183e-03*pow(log(x),4)"},}},
  {"hadHcalp3_nhf",          {{"ispositive", "0"}, {"initial", "0"}, {"freeze", "0"}, {"type","hadHcalp3"},          {"appliesTo","nhf"},  {"form", "-4.683e+00+3.771e+00*log(x)-1.038e+00*pow(log(x),2)+1.224e-01*pow(log(x),3)-5.071e-03*pow(log(x),4)"},}},
  {"hadHcalp3_nef",          {{"ispositive", "0"}, {"initial", "0"}, {"freeze", "0"}, {"type","hadHcalp3"},          {"appliesTo","nef"},  {"form", "+8.141e+00-6.447e+00*log(x)+1.779e+00*pow(log(x),2)-2.020e-01*pow(log(x),3)+7.952e-03*pow(log(x),4)"},}},

  {"hadHcalZB097",           {{"ispositive", "0"}, {"initial", "0"}, {"freeze", "0"}, {"type","hadHcalZB097"},       {"appliesTo","Resp"}, {"form", "-7.714e+00+7.581e+00*log(x)-2.665e+00*pow(log(x),2)+3.940e-01*pow(log(x),3)-2.049e-02*pow(log(x),4)"},}},
  {"hadHcalZB097_chf",       {{"ispositive", "0"}, {"initial", "0"}, {"freeze", "0"}, {"type","hadHcalZB097"},       {"appliesTo","chf"},  {"form", "+3.089e+00-2.918e+00*log(x)+1.006e+00*pow(log(x),2)-1.477e-01*pow(log(x),3)+7.648e-03*pow(log(x),4)"},}},
  {"hadHcalZB097_nhf",       {{"ispositive", "0"}, {"initial", "0"}, {"freeze", "0"}, {"type","hadHcalZB097"},       {"appliesTo","nhf"},  {"form", "-9.328e+00+7.999e+00*log(x)-2.493e+00*pow(log(x),2)+3.336e-01*pow(log(x),3)-1.605e-02*pow(log(x),4)"},}},
  {"hadHcalZB097_nef",       {{"ispositive", "0"}, {"initial", "0"}, {"freeze", "0"}, {"type","hadHcalZB097"},       {"appliesTo","nef"},  {"form", "+4.981e+00-4.021e+00*log(x)+1.166e+00*pow(log(x),2)-1.446e-01*pow(log(x),3)+6.489e-03*pow(log(x),4)"},}},

  {"hadHcalZB100",           {{"ispositive", "0"}, {"initial", "0"}, {"freeze", "0"}, {"type","hadHcalZB100"},       {"appliesTo","Resp"}, {"form", "-9.705e+00+8.983e+00*log(x)-2.986e+00*pow(log(x),2)+4.190e-01*pow(log(x),3)-2.066e-02*pow(log(x),4)"},}},
  {"hadHcalZB100_chf",       {{"ispositive", "0"}, {"initial", "0"}, {"freeze", "0"}, {"type","hadHcalZB100"},       {"appliesTo","chf"},  {"form", "+4.857e+00-4.366e+00*log(x)+1.427e+00*pow(log(x),2)-1.990e-01*pow(log(x),3)+9.812e-03*pow(log(x),4)"},}},
  {"hadHcalZB100_nhf",       {{"ispositive", "0"}, {"initial", "0"}, {"freeze", "0"}, {"type","hadHcalZB100"},       {"appliesTo","nhf"},  {"form", "-9.351e+00+7.887e+00*log(x)-2.408e+00*pow(log(x),2)+3.141e-01*pow(log(x),3)-1.461e-02*pow(log(x),4)"},}},
  {"hadHcalZB100_nef",       {{"ispositive", "0"}, {"initial", "0"}, {"freeze", "0"}, {"type","hadHcalZB100"},       {"appliesTo","nef"},  {"form", "+3.789e+00-2.934e+00*log(x)+8.056e-01*pow(log(x),2)-9.292e-02*pow(log(x),3)+3.776e-03*pow(log(x),4)"},}},

  {"hadHcalZB106",           {{"ispositive", "0"}, {"initial", "0"}, {"freeze", "0"}, {"type","hadHcalZB106"},       {"appliesTo","Resp"}, {"form", "-2.023e+01+1.720e+01*log(x)-5.226e+00*pow(log(x),2)+6.681e-01*pow(log(x),3)-2.976e-02*pow(log(x),4)"},}},
  {"hadHcalZB106_chf",       {{"ispositive", "0"}, {"initial", "0"}, {"freeze", "0"}, {"type","hadHcalZB106"},       {"appliesTo","chf"},  {"form", "+1.035e+01-8.927e+00*log(x)+2.779e+00*pow(log(x),2)-3.677e-01*pow(log(x),3)+1.719e-02*pow(log(x),4)"},}},
  {"hadHcalZB106_nhf",       {{"ispositive", "0"}, {"initial", "0"}, {"freeze", "0"}, {"type","hadHcalZB106"},       {"appliesTo","nhf"},  {"form", "-1.206e+01+9.880e+00*log(x)-2.898e+00*pow(log(x),2)+3.578e-01*pow(log(x),3)-1.535e-02*pow(log(x),4)"},}},
  {"hadHcalZB106_nef",       {{"ispositive", "0"}, {"initial", "0"}, {"freeze", "0"}, {"type","hadHcalZB106"},       {"appliesTo","nef"},  {"form", "+3.033e+00-2.060e+00*log(x)+4.513e-01*pow(log(x),2)-3.274e-02*pow(log(x),3)+1.314e-04*pow(log(x),4)"},}},

  {"ecalm3",                 {{"ispositive", "0"}, {"initial", "0"}, {"freeze", "0"}, {"type","ecalm3"},             {"appliesTo","Resp"}, {"form", "-9.628e-01+1.009e+00*log(x)-4.167e-01*pow(log(x),2)+5.326e-02*pow(log(x),3)-2.206e-03*pow(log(x),4)"},}},
  {"ecalm3_chf",             {{"ispositive", "0"}, {"initial", "0"}, {"freeze", "0"}, {"type","ecalm3"},             {"appliesTo","chf"},  {"form", "+5.813e+00-4.086e+00*log(x)+1.102e+00*pow(log(x),2)-1.220e-01*pow(log(x),3)+4.691e-03*pow(log(x),4)"},}},
  {"ecalm3_nhf",             {{"ispositive", "0"}, {"initial", "0"}, {"freeze", "0"}, {"type","ecalm3"},             {"appliesTo","nhf"},  {"form", "-1.636e-01+1.183e-01*log(x)-9.015e-03*pow(log(x),2)-3.084e-03*pow(log(x),3)+3.948e-04*pow(log(x),4)"},}},
  {"ecalm3_nef",             {{"ispositive", "0"}, {"initial", "0"}, {"freeze", "0"}, {"type","ecalm3"},             {"appliesTo","nef"},  {"form", "-5.823e+00+4.103e+00*log(x)-1.133e+00*pow(log(x),2)+1.300e-01*pow(log(x),3)-5.306e-03*pow(log(x),4)"},}},

  {"ecalGain1p3",            {{"ispositive", "0"}, {"initial", "0"}, {"freeze", "0"}, {"type","ecalGain1p3"},        {"appliesTo","Resp"}, {"form", "-8.351e+00+6.497e+00*log(x)-1.773e+00*pow(log(x),2)+1.998e-01*pow(log(x),3)-7.719e-03*pow(log(x),4)"},}},
  {"ecalGain1p3_chf",        {{"ispositive", "0"}, {"initial", "0"}, {"freeze", "0"}, {"type","ecalGain1p3"},        {"appliesTo","chf"},  {"form", "+2.612e+00-2.118e+00*log(x)+6.084e-01*pow(log(x),2)-7.302e-02*pow(log(x),3)+3.065e-03*pow(log(x),4)"},}},
  {"ecalGain1p3_nhf",        {{"ispositive", "0"}, {"initial", "0"}, {"freeze", "0"}, {"type","ecalGain1p3"},        {"appliesTo","nhf"},  {"form", "-9.112e-01+8.269e-01*log(x)-2.729e-01*pow(log(x),2)+3.890e-02*pow(log(x),3)-2.026e-03*pow(log(x),4)"},}},
  {"ecalGain1p3_nef",        {{"ispositive", "0"}, {"initial", "0"}, {"freeze", "0"}, {"type","ecalGain1p3"},        {"appliesTo","nef"},  {"form", "-3.583e-01+2.078e-01*log(x)-2.141e-02*pow(log(x),2)-4.756e-03*pow(log(x),3)+7.013e-04*pow(log(x),4)"},}},

  {"ecalGain6p3",            {{"ispositive", "0"}, {"initial", "0"}, {"freeze", "0"}, {"type","ecalGain6p3"},        {"appliesTo","Resp"}, {"form", "-5.587e+00+4.992e+00*log(x)-1.596e+00*pow(log(x),2)+2.157e-01*pow(log(x),3)-1.033e-02*pow(log(x),4)"},}},
  {"ecalGain6p3_chf",        {{"ispositive", "0"}, {"initial", "0"}, {"freeze", "0"}, {"type","ecalGain6p3"},        {"appliesTo","chf"},  {"form", "+2.414e+00-2.133e+00*log(x)+6.755e-01*pow(log(x),2)-9.053e-02*pow(log(x),3)+4.308e-03*pow(log(x),4)"},}},
  {"ecalGain6p3_nhf",        {{"ispositive", "0"}, {"initial", "0"}, {"freeze", "0"}, {"type","ecalGain6p3"},        {"appliesTo","nhf"},  {"form", "+7.040e-01-5.838e-01*log(x)+1.728e-01*pow(log(x),2)-2.155e-02*pow(log(x),3)+9.501e-04*pow(log(x),4)"},}},
  {"ecalGain6p3_nef",        {{"ispositive", "0"}, {"initial", "0"}, {"freeze", "0"}, {"type","ecalGain6p3"},        {"appliesTo","nef"},  {"form", "-3.593e+00+3.112e+00*log(x)-9.667e-01*pow(log(x),2)+1.272e-01*pow(log(x),3)-5.959e-03*pow(log(x),4)"},}},

  {"ecalGain12p3",           {{"ispositive", "0"}, {"initial", "0"}, {"freeze", "0"}, {"type","ecalGain12p3"},       {"appliesTo","Resp"}, {"form", "+1.542e+01-1.289e+01*log(x)+3.895e+00*pow(log(x),2)-4.819e-01*pow(log(x),3)+2.084e-02*pow(log(x),4)"},}},
  {"ecalGain12p3_chf",       {{"ispositive", "0"}, {"initial", "0"}, {"freeze", "0"}, {"type","ecalGain12p3"},       {"appliesTo","chf"},  {"form", "-1.059e+01+8.124e+00*log(x)-2.321e+00*pow(log(x),2)+2.773e-01*pow(log(x),3)-1.169e-02*pow(log(x),4)"},}},
  {"ecalGain12p3_nhf",       {{"ispositive", "0"}, {"initial", "0"}, {"freeze", "0"}, {"type","ecalGain12p3"},       {"appliesTo","nhf"},  {"form", "-9.448e-02+1.776e-02*log(x)-1.417e-03*pow(log(x),2)-4.896e-04*pow(log(x),3)+5.913e-05*pow(log(x),4)"},}},
  {"ecalGain12p3_nef",       {{"ispositive", "0"}, {"initial", "0"}, {"freeze", "0"}, {"type","ecalGain12p3"},       {"appliesTo","nef"},  {"form", "+1.063e+01-8.107e+00*log(x)+2.314e+00*pow(log(x),2)-2.758e-01*pow(log(x),3)+1.158e-02*pow(log(x),4)"},}},

  {"trkEff0999Nm1",          {{"ispositive", "0"}, {"initial", "0"}, {"freeze", "0"}, {"type","trkEff0999Nm1"},      {"appliesTo","Resp"}, {"form", "-2.415e+00+2.029e+00*log(x)-6.003e-01*pow(log(x),2)+7.296e-02*pow(log(x),3)-3.117e-03*pow(log(x),4)"},}},
  {"trkEff0999Nm1_chf",      {{"ispositive", "0"}, {"initial", "0"}, {"freeze", "0"}, {"type","trkEff0999Nm1"},      {"appliesTo","chf"},  {"form", "-1.755e+00+9.661e-01*log(x)-8.482e-02*pow(log(x),2)-2.171e-02*pow(log(x),3)+2.382e-03*pow(log(x),4)"},}},
  {"trkEff0999Nm1_nhf",      {{"ispositive", "0"}, {"initial", "0"}, {"freeze", "0"}, {"type","trkEff0999Nm1"},      {"appliesTo","nhf"},  {"form", "-2.462e+00+2.152e+00*log(x)-6.830e-01*pow(log(x),2)+9.276e-02*pow(log(x),3)-4.384e-03*pow(log(x),4)"},}},
  {"trkEff0999Nm1_nef",      {{"ispositive", "0"}, {"initial", "0"}, {"freeze", "0"}, {"type","trkEff0999Nm1"},      {"appliesTo","nef"},  {"form", "+5.553e+00-4.180e+00*log(x)+1.069e+00*pow(log(x),2)-1.073e-01*pow(log(x),3)+3.559e-03*pow(log(x),4)"},}},

  {"trkEff0998Nm1",          {{"ispositive", "0"}, {"initial", "0"}, {"freeze", "0"}, {"type","trkEff0998Nm1"},      {"appliesTo","Resp"}, {"form", "-4.645e+00+3.917e+00*log(x)-1.161e+00*pow(log(x),2)+1.410e-01*pow(log(x),3)-6.012e-03*pow(log(x),4)"},}},
  {"trkEff0998Nm1_chf",      {{"ispositive", "0"}, {"initial", "0"}, {"freeze", "0"}, {"type","trkEff0998Nm1"},      {"appliesTo","chf"},  {"form", "-5.411e+00+3.418e+00*log(x)-5.853e-01*pow(log(x),2)+6.119e-03*pow(log(x),3)+2.638e-03*pow(log(x),4)"},}},
  {"trkEff0998Nm1_nhf",      {{"ispositive", "0"}, {"initial", "0"}, {"freeze", "0"}, {"type","trkEff0998Nm1"},      {"appliesTo","nhf"},  {"form", "-5.007e+00+4.371e+00*log(x)-1.388e+00*pow(log(x),2)+1.887e-01*pow(log(x),3)-8.947e-03*pow(log(x),4)"},}},
  {"trkEff0998Nm1_nef",      {{"ispositive", "0"}, {"initial", "0"}, {"freeze", "0"}, {"type","trkEff0998Nm1"},      {"appliesTo","nef"},  {"form", "+1.121e+01-8.435e+00*log(x)+2.159e+00*pow(log(x),2)-2.174e-01*pow(log(x),3)+7.266e-03*pow(log(x),4)"},}},

  {"trkEffNtrk1m3",          {{"ispositive", "0"}, {"initial", "0"}, {"freeze", "0"}, {"type","trkEffNtrk1m3"},      {"appliesTo","Resp"}, {"form", "+6.020e+00-5.647e+00*log(x)+1.649e+00*pow(log(x),2)-1.955e-01*pow(log(x),3)+8.229e-03*pow(log(x),4)"},}},
  {"trkEffNtrk1m3_chf",      {{"ispositive", "0"}, {"initial", "0"}, {"freeze", "0"}, {"type","trkEffNtrk1m3"},      {"appliesTo","chf"},  {"form", "+4.784e+00-5.104e+00*log(x)+1.600e+00*pow(log(x),2)-1.994e-01*pow(log(x),3)+8.710e-03*pow(log(x),4)"},}},
  {"trkEffNtrk1m3_nhf",      {{"ispositive", "0"}, {"initial", "0"}, {"freeze", "0"}, {"type","trkEffNtrk1m3"},      {"appliesTo","nhf"},  {"form", "-2.550e+00+2.776e+00*log(x)-8.980e-01*pow(log(x),2)+1.156e-01*pow(log(x),3)-5.216e-03*pow(log(x),4)"},}},
  {"trkEffNtrk1m3_nef",      {{"ispositive", "0"}, {"initial", "0"}, {"freeze", "0"}, {"type","trkEffNtrk1m3"},      {"appliesTo","nef"},  {"form", "-3.124e+00+3.086e+00*log(x)-9.373e-01*pow(log(x),2)+1.149e-01*pow(log(x),3)-4.982e-03*pow(log(x),4)"},}},

  {"trkEffNtrk2ToInfm3",     {{"ispositive", "0"}, {"initial", "0"}, {"freeze", "0"}, {"type","trkEffNtrk2ToInfm3"}, {"appliesTo","Resp"}, {"form", "+3.457e+00-1.972e+00*log(x)+2.966e-01*pow(log(x),2)-9.399e-03*pow(log(x),3)-4.569e-04*pow(log(x),4)"},}},
  {"trkEffNtrk2ToInfm3_chf", {{"ispositive", "0"}, {"initial", "0"}, {"freeze", "0"}, {"type","trkEffNtrk2ToInfm3"}, {"appliesTo","chf"},  {"form", "-8.492e+00+7.806e+00*log(x)-2.447e+00*pow(log(x),2)+2.963e-01*pow(log(x),3)-1.222e-02*pow(log(x),4)"},}},
  {"trkEffNtrk2ToInfm3_nhf", {{"ispositive", "0"}, {"initial", "0"}, {"freeze", "0"}, {"type","trkEffNtrk2ToInfm3"}, {"appliesTo","nhf"},  {"form", "-4.057e+00+2.813e+00*log(x)-6.923e-01*pow(log(x),2)+7.966e-02*pow(log(x),3)-3.445e-03*pow(log(x),4)"},}},
  {"trkEffNtrk2ToInfm3_nef", {{"ispositive", "0"}, {"initial", "0"}, {"freeze", "0"}, {"type","trkEffNtrk2ToInfm3"}, {"appliesTo","nef"},  {"form", "+1.414e+01-1.197e+01*log(x)+3.551e+00*pow(log(x),2)-4.294e-01*pow(log(x),3)+1.815e-02*pow(log(x),4)"},}},

};
