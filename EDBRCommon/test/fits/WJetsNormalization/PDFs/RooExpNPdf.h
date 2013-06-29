/*****************************************************************************
 * Project: RooFit                                                           *
 *                                                                           *
  * This code was autogenerated by RooClassFactory                            * 
 *****************************************************************************/

#ifndef ROOEXPNPDF
#define ROOEXPNPDF

#include "RooAbsPdf.h"
#include "RooRealProxy.h"
#include "RooCategoryProxy.h"
#include "RooAbsReal.h"
#include "RooAbsCategory.h"
 
class RooExpNPdf : public RooAbsPdf {
public:
  RooExpNPdf() {} ; 
  RooExpNPdf(const char *name, const char *title,
	      RooAbsReal& _x,
	      RooAbsReal& _c,
	      RooAbsReal& _n);
  RooExpNPdf(const RooExpNPdf& other, const char* name=0) ;
  virtual TObject* clone(const char* newname) const { return new RooExpNPdf(*this,newname); }
  inline virtual ~RooExpNPdf() { }

protected:

  RooRealProxy x ;
  RooRealProxy c ;
  RooRealProxy n ;
  
  Double_t evaluate() const ;

private:

  ClassDef(RooExpNPdf,1) // Your description goes here...
};
 
#endif