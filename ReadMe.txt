A linear multiple actuator cylinder code, an initial release of the planned full multiple actuator cylinder code.

Produced by Angus Buck, PhD student at the University of Strathclyde, following the work of Andrew Ning (aning@byu.edu) in "Actuator cylinder theory for multiple vertical axis wind turbines" (doi:10.5194/wes-1-327-2016).

This version (V0.1.1) is only linear and is known to not be fully accurate; it is however certainly giving results in the right ballpark, and will continue to be updated as I iron out these issues. V0.1.2 will be the same code with an added folder showing validations, and the upcoming updates are expected to be inclusion of nonlinearities and an improvement of 'fsolve' converging.
Email angus.buck.2016@uni.strath.ac.uk if you want to be added to a mailing list for when this code is updated :-)

Openess and honestly: If you know how to work Julia and are comfortable using a less-commented code, it's very likely that using Andrew Ning's own code, found at https://github.com/byuflowlab/vawt-ac, will give more accurate results than this code in its current state.
If however you want a more "pick-up-and-play" model, this code is very thoroughly commented and easy to use: just open the file, set your variables, change your aerofoil data if you want to, and hit "run". This code also has batch-running capability and the ability to compare results with different inputs, which I'm not sure if Ning's has (I don't speak Julia!)

I believe there are other good actuator cylinder codes (references below), but aside from HAWC2, I'm not as sure on their availability.
Edgar Martinez-Ojeda et al: doi.org/10.5194/wes-6-1061-2021
HAWC2: Madsen et al, Implementation of the Actuator Cylinder Flow Model in the HAWC2 code for Aeroelastic Simulations on Vertical Axis Wind Turbines
Ang Li: Double actuator cylinder (AC) model of a tandem vertical-axis wind turbine (VAWT) counter-rotating rotor concept operating in different wind conditions (2017)
Zhengshun Cheng: doi.org/10.1016/j.egypro.2016.09.232. 