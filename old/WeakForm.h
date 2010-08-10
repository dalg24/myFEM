#ifndef WEAKFORM_H
#define WEAKFORM_H

#include <iostream>

/** 
* @class WeakForm
*
* @brief This is a class for the weak form.
*
*
* @author DaLG <lebrungd@neo.tamu.edu>
*/

class WeakForm {

public:
  /// Constructor.
  WeakForm();
 
  /// Copy constructor.
  WeakForm(const WeakForm&);

  /// Destructor.  
  ~WeakForm();

  WeakForm& operator=(const WeakForm&);

  ///
  void set_bilinear_form(const BilinearForm) const;

  ///
  void set_linear_form(const LinearForm) const;

private:
  /// Bilinear form.
  BilinearForm *bilinear_form;

  /// Linear form.
  LinearForm *linear_form;

};

#endif // WEAKFORM_H
