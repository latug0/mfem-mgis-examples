# Material Models and Constitutive Laws

This file contains the constitutive laws implemented via **TFEL/MFront** for the thermomechanical coupling simulation. The system models two distinct materials: the cladding/structure (**ALFENI**) and the fuel (**U3Si2**).

---

## Material 1: ALFENI (Cladding / Structure)

ALFENI is modeled with an elasto-plastic behavior coupled with temperature, accounting for large strains and solving the heat equation after considering the effects of
the deformation.

### Mechanical Behavior: Plasticity with Hardening (`_IsotropicLinearHardeningPlasticity`)
This behavior implements Von Mises plasticity with linear isotropic hardening. It is formulated in large strains using the **Hencky** strain measure.

* **Thermal Expansion:** The thermal expansion coefficient is set to $15 \times 10^{-6}\text{ K}^{-1}$ with a reference temperature of $293.15\text{ K}$.
* **Exported Variables (Post-processing):** The temperature, yield strength ($\sigma_0$), and hardening modulus ($H$) are exported as auxiliary state variables to facilitate visualization (e.g., in ParaView).

### Material Laws (ALFENI)
The mechanical properties of ALFENI depend on temperature according to the following laws:

* **Young's Modulus & Poisson's Ratio:** Both the Young's modulus and Poisson's ratio are assumed to be **constant**.
* **Yield Strength (`YieldStrength`):**
  * Constant at $200\text{ MPa}$ for $T < 573.15\text{ K}$.
  * Linear thermal softening above $573.15\text{ K}$ (slope of $-200\text{ kPa/K}$).
  * Floor value set at $10\text{ MPa}$.
* **Hardening Modulus (`HardeningModulus`):** (Defined via the C++ header `ALFENI_HardeningModulus.hxx`)
  * Constant at $1500\text{ MPa}$ up to $573.15\text{ K}$.
  * Linear thermal softening above this threshold (slope of $-1.5\text{ MPa/K}$).
  * Floor value set at $100\text{ MPa}$.

### Coupled Thermal Behavior (`_ThermiqueCouplee`)
The heat equation is solved taking into account the geometry deformation.
* The effective conductivity $K_{eff}$ depends on the deformation gradient $F$ (via the right Cauchy-Green tensor $C$ and the Jacobian $J$).
* The intrinsic conductivity $k_0$ evolves linearly with temperature: $k_0 = A + B(T - 273.15)$, with $A = 121.0$ and $B = 0.15$.
* The exact tangent operator (both thermal and mechanical parts) is provided to the solver to ensure optimal convergence.

---

## Material 2: U3Si2 (Fuel)

The U3Si2 fuel is subjected to in-pile irradiation phenomena. Its mechanical behavior is dominated by irradiation creep and solid swelling induced by fissions.

### Solid Swelling Model (`_SolidSwelling`)
This point-wise model computes the evolution of the solid volumetric swelling $S$ based on the local power density.

* **Inputs:** External power density (`Pow`) and average energy released per fission ($fe = 200\text{ MeV}$).
* **Integration:** The swelling rate is integrated implicitly using the average power density over the time step.
* **Rate:** Proportional constant set to $6.2 \times 10^{-29}$ (adjusted via energy unit conversions).

### Mechanical Behavior: Irradiation Creep (`_NortonPRQ`)
Unlike ALFENI, U3Si2 does not undergo classical plasticity but instead creeps under irradiation (modified Norton-type law). This model also uses the **Hencky** strain measure for large deformations.

* **Coupling with Swelling:** The swelling $sw$ computed by the `_SolidSwelling` model is passed to this mechanical law as an external state variable. It is converted into an isotropic inelastic strain tensor.
* **Flow Rule (Creep):** 
  * Driven by the fission rate ($fr = Pow / E_f$), where $E_f = 3.204 \times 10^{-11}\text{ J/fission}$.
  * Creep is proportional to the Von Mises equivalent stress ($seq$) multiplied by the fission rate and a constant $A = 500 \times 10^{-37}$.
* **Post-processing:** The scalar swelling is assigned to an auxiliary variable (`SwellingExport` of type `strain`) to force its export to visualization tools.