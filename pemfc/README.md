# PEMFC Model

This model simulates constant-current operation of a PEM fuel cell.  At present, the domain includes:
- Anode and cathode gas diffusion layer (GDL), with constant gas-phase composition (`C_k`) and electric potential(`phi`)
- Anode and cathode catalyst layer (CL), which includes carbon-supported Pt catalyst phase, coated by Nafion ionomer (which provides an ion conduction pathway to the membrane), and gas phase.

The gas-phase composition and double-layer electric potential (`dPhi_dl`) in the catalyst layers are simulated via physically-based partial differential equations.  Modeled phenomena include:
- Global half-cell electrochemical reactions, with associated Faradaic current (`i_Far`).
- Double layer charging current (`i_dl`) to maintain charge neutrality.
- Ionic (`i_io`) and electronic (`i_el`) current into and out of the CL, set equal to `i_ext`.
- Gas-phase transport (diffusion + convection) between the GDL and CL.

The model can be run in one of three modes:
1. Call the model directly for a single external current density and temperature, by running `pemfc_model.py`
2. Run a polarization curve (which calls the model for a range of current densities), by running `pemfc_polarization.py`
3. Run polarization curves for a range of temperatures, by running `pemfc_temp_study.py`.  This calls the polarization curve function for a range of temperatures.

Setting parameters and conditions:
The full range of parmeters and input values can be set in `pemfc_inputs.py`. You can overwrite two values from this file, when running the model - the external current and the temperature.  Doing this is demonstrated by the function `pemfc_polarization.py`.
