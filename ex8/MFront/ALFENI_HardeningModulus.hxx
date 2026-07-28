#pragma once 

inline double computeHardeningModulus(const double T_actuelle) {
  double H = 1500e6; 
  
  if (T_actuelle > 573.15) {
    H -= 1.5e6 * (T_actuelle - 573.15);
  }
  
  if (H < 100e6) {
    H = 100e6; 
  }
  
  return H;
}