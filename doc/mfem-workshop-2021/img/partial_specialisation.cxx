template<typename T,typename T2>
typename std::enable_if<
  ((tfel::meta::Implements<T,StensorConcept>::cond)  && (StensorTraits<T>::dime==3u)&&
   (tfel::meta::Implements<T2,StensorConcept>::cond) && (StensorTraits<T2>::dime==3u)&&
   (tfel::typetraits::IsFundamentalNumericType<StensorNumType<T2>>::cond)),
  stensor<3u,StensorNumType<T>>
 >::type
convertCorotationnalCauchyStressToSecondPiolaKirchhoffStress(const T& s, const T2& U){
  using real = tfel::typetraits::base_type<StensorNumType<T>>;
  constexpr real cste = Cste<real>::sqrt2;
  const auto J  = det(U);
  const auto iU = invert(U);
  return {J*(s[2]*iU[4]*iU[4]+(cste*s[5]*iU[3]+2*s[4]*iU[0])*iU[4]+s[1]*iU[3]*iU[3]+2*s[3]*iU[0]*iU[3]+2*s[0]*iU[0]*iU[0])/2,
	  J*(s[2]*iU[5]*iU[5]+(cste*s[4]*iU[3]+2*s[5]*iU[1])*iU[5]+s[0]*iU[3]*iU[3]+2*s[3]*iU[1]*iU[3]+2*s[1]*iU[1]*iU[1])/2,
	  J*(s[1]*iU[5]*iU[5]+(cste*s[3]*iU[4]+2*s[5]*iU[2])*iU[5]+s[0]*iU[4]*iU[4]+2*s[4]*iU[2]*iU[4]+2*s[2]*iU[2]*iU[2])/2,
	  J*((cste*s[2]*iU[4]+s[5]*iU[3]+cste*s[4]*iU[0])*iU[5]+(s[4]*iU[3]+cste*s[5]*iU[1])*iU[4]+s[3]*iU[3]*iU[3]+(2*s[1]*iU[1]+2*s[0]*iU[0])*iU[3]+2*s[3]*iU[0]*iU[1])/2,
	  J*((s[5]*iU[4]+cste*s[1]*iU[3]+cste*s[3]*iU[0])*iU[5]+s[4]*iU[4]*iU[4]+(s[3]*iU[3]+2*s[2]*iU[2]+2*s[0]*iU[0])*iU[4]+cste*s[5]*iU[2]*iU[3]+2*s[4]*iU[0]*iU[2])/2,
	  J*(s[5]*iU[5]*iU[5]+(s[4]*iU[4]+s[3]*iU[3]+2*s[2]*iU[2]+2*s[1]*iU[1])*iU[5]+(cste*s[0]*iU[3]+cste*s[3]*iU[1])*iU[4]+cste*s[4]*iU[2]*iU[3]+2*s[5]*iU[1]*iU[2])/2};
}
