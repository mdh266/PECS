*VIEW ON https://mdh266.github.io/PECS/*

\mainpage  

\author Michael Harmon

\section Overview

\subsection Introduction

Large-scale utilization of photovoltaic (PV) devices, or solar cells, has been hampered 
for years due to high costs and lack of energy storage mechanisms.  
<a href="https://en.wikipedia.org/wiki/Photoelectrochemical_cell">Photoelectrochemical solar cells</a>  
(PECs), are an attractive alternative to conventional solid state PV devices.  
PECs such as those depicted below 
are able to directly convert solar energy into hydrogen fuel.  The hydrogen fuel can then 
be used at a later time to generate electricity. The typical setup of a PEC is shown below

\image html PEC.png

A PEC consists of four main components: the solid semiconductor electrode, 
the liquid electrolyte component, the semiconductor-electrolyte interface and the counter 
(metal or semiconductor) electrode.  When sunlight shines on the semiconductor component, 
photons are absorbed and generate electron-hole pairs. These electrons and holes are 
separated by a built-in electric field within the semiconductor. The separation of 
the electrons and holes leads to an electrical current in the cell and the accumulation 
of charges at the semiconductor-electrolyte interface.  At the interface, the photo-generated 
electrons or holes induce a chemical reaction at the semiconductor electrode.
A similar chemical reaction also occurs at the counter electrode. These chemical reactions 
create or eliminate reductant-oxidant (redox) species in the electrolyte leading
to the generation of hydrogen fuel.  

Research on PECs has traditionally focused on planar cell designs, but recently there has 
been interest in cell designs that use thin nanostructured wires such as those depicted below.

\image html planar_vs_wire-eps-converted-to.png

In a planar device such, photo-generated electrons and holes are collected
in directions parallel to photon absorption. In order for PECs' to achieve sufficient 
energy conversion efficiencies to be commercially viable the electron/hole diffusion length 
(\f$L_{D}\f$) (the average distance an electron/hole can be travel without being eliminated)
must be larger than the absorption length (\f$1/\alpha\f$) (the average distance a photon
will penetrate the semiconductor crystal before generating an electron-hole pair). 
This constraint necessitates the use of expensive, high quality crystals that either have 
large diffusion lengths or small absorption lengths. In PECs that use a nanowire design the electron/hole collection and photon absorption directions are 
decoupled, thereby alleviating the need for high quality crystals to attain sufficient energy 
conversion efficiencies.  Scientists believe that the shape of these nanowires can effect 
the energy conversion efficiencies (see <a href="http://pubs.rsc.org/en/content/articlelanding/2012/ee/c1ee02518j#!divAbstract">here</a> for more information).
This software therefore allows nanowires which are more cylindrical or conic like those seen below:


\image html Cylinder.png



This software is designed to simulate the dynamics of the reactive interface between a
semiconductor and electrolyte in a 2D. The nanowire solar cell is then considered to be axially symmetric in order
to reduce the 3D nanowire to a 2D domain. The interface of between the semiconductor and electrolyte
make up a "half cell" of a 
<a href="https://en.wikipedia.org/wiki/Photoelectrochemical_cell">photoelectrochemical cell</a>.
The main challenges in constructing numerical algorithms that produces reliable simulations of PECs
are due to the highly nonlinear nature of the semiconductor and electrolyte systems as well 
as the nonlinear coupling between the two systems at the interface. In addition, the 
evolution problem under consideration is effectively multi-scale in the sense that the evolution 
of the system in the semiconductor and the corresponding one in the electrolyte evolve at 
different time scales due to the quantitative scaling differences in their relevant physical parameters. 
Furthermore, regions of stiffness caused by boundary layer formation where sharp transitions in  
densities and electric potential occur near the interface and pose severe constraints on the choice of 
discretization strategy in order to maintain numerical stability.


\subsection Requirements 
The requirements for this software are 
<a href="https://www.dealii.org">deal.ii library</a> version 8.3.0 or higher,
and <a href="https://www.cmake.org">CMake</a> version 2.8 or higher.  The code will
automatically run in parallel using the <a href="https://www.threadingbuildingblocks.org/">
Thread Building Blocks</a>.
See
<a href="https://www.dealii.org/developer/doxygen/deal.II/group__threads.html">
deal.ii's explanation on parallel computing with shared memory</a> to see why this 
necessary and how it works. 


\subsection Using the code
The source code can be downloaded from <a href="https://github.com/mdh266/PECS">here</a>.
You first need to obtain and install a copy of the 
<a href="https://www.dealii.org">deal.ii library</a> 
version 8.3.0 or higher. After downloading and installing the deal.II library.  
cd into the PECS directory.  

To generate a make file to compile the source code run the following command in the <code>PECS/</code> directory,

<code> cmake . -DDEAL_II_DIR="path to deal.II library" </code>

On a mac, if you downloaded the binaries of deal.ii library instead run:

<code> cmake . </code>

Once the CMake build is complete you can run the command:

<code>make release</code>

to compile the code.

To run the code use the command:

<code>./main</code>

The resulting output files of the simulation are in <a href="http://www.vtk.org/"> VTK </a> format and can 
viewed using <a href="http://www.paraview.org/"> Paraview</a>.

\subsection Parameters

Material and design choices can be chosen by the user through the input file
<code>input_file.prm</code>.

Most of the parameters in <code>input_file.prm</code> will be self explanatory after
reading the Model section and from the comments in the file. However, we explain some of them
here for clarity.  There are two end times in <code>input_file.prm</code> in case in the initial
simulation has not fully converged to steady state. Instead of restarting the simulation from the 
initial conditions and rerunning the simulation we can restart the simulation and use the 
end conditions as our new starting conditions; this is achieved by setting:

<code> set restart status     = true</code>

And setting the <code>end time 2</code> to the desired new end time.

A typical mesh for the semiconductor is seen below,

\image html semiconductor-grid.png

The height of the mesh is the height of the nanowire and is set in,

<code>  set mesh height    = ... </code>

The top width is set by setting the top radius of the nano wire,

<code>  set radius one     =  ... </code>

the bottom width is set by setting the bottom radius of the nano wire,

<code>  set radius two     =  ... </code>

Refining the mesh globally is achieved by,

<code> set global refinements = ... </code>

Refining the mesh locally is achieved by,

<code> set local refinements  = ... </code>

Local refinement only occurs within a distance of the interface (on the right) and this distance is set by,

<code> set boundary layer = ... </code>

The whole simulation mesh is seen below, 

\image html Poisson-grid.png

The details of this mesh (refinements, height, radii, boundary layer) are controlled as described above. 
The entire width of the domain is set by,

<code>  set mesh length    = ... </code>

Both radii **MUST** be smaller than the mesh width. See Grid_Maker::Grid for more information.  
Finally, the left and the right boundaries of the the entire domain will have Dirchlet
boundary conditions (explained in Model section), while the top and bottom of the domain
will have Neumann conditions (explained in Model section) by,


<code>  set insulated    = true </code>

otherwise the top and bottom of the domain will have Dirichlet conditions.  Finally, the
initial conditions are set in the source code file:

<code> InitialConditions.cpp </code>

and must be set before compilation.

**Note**

If your simulation seems to be losing stability, there are two options:

1.) Choose a smaller time step, in 

<code>  set time step size     = ... </code>

2.) Refining the mesh.



\section Background

The software provides simulations of a photoelectrocrhemical solar cell
by solving the reaction diffusion equations that describe the macroscopic dynamics of 
charge transport in photoelectrochemical (PEC) solar cell.  
The main objective is to accurately capture the reactive dynamics of the 
semiconductor-electrolyte interface. The underlying spatial descritizations
are based on the local discontinuous Galerkin (LDG_System::LDG) 
and mixed finite element method (MixedPoisson::MixedFEM).  The use of specific tailored
implicit-explicit (IMEX) methods to capture transient solutions of the
time dependent nonlinear systems of PDES. The main class in which everything is put 
together is in SOLARCELL::SolarCellProblem. Information on the numerical methods are
provided in the links above as well as in the 
<a href="http://www.ma.utexas.edu/users/gamba/papers/Harmon-IMG-Ren.pdf"> paper </a>
based of this software.  


We present the mathematical model below, however a more through introduction can be 
found <a href="http://www.ma.utexas.edu/users/gamba/papers/Semi-Elec-Interf.pdf"> here 
</a>.


\subsection Model
We focus on the reactive dynamics of the semiconductor-electrolyte interface.
Therfore our domain is the half cell, whose abstract representation is, 
\image html DomainDecomposition.png
	

In the semiconductor component (\f$ \Omega_{S} \f$) the transport of electrons 
\f$(\rho_{n})\f$ and holes \f$(\rho_{p})\f$ is governed by 
the drift-diffusion-Poisson system of equations for electrons
\f[ \begin{align} 
\frac{\partial \rho_{n}}{\partial t}
\, + \, 
\boldsymbol \nabla  \cdot  \,\left( - \alpha_{n} \, \mu_{n} \boldsymbol \nabla \Phi \, \rho_{n} 
\, - \,
D_{n} \  \boldsymbol \nabla \, \rho_{n}  \, \right) 
\; &= \; 
R(\rho_{n}, \rho_{p})
\, + \, 
G && \text{in} \ (0,T] \ \times \ \Omega_{S} \nonumber \\
\frac{\partial \rho_{p}}{\partial t} 
\, + \, 
\boldsymbol \nabla \cdot  \, \left( - \alpha_{p} \, \mu_{p} \boldsymbol \nabla \Phi \, \rho_{p} 
\, - \, 
D_{p} \  \boldsymbol \nabla \, \rho_{p} \, \right) 
\; &= \;
R(\rho_{n}, \rho_{p})
\, + \, 
G
&& \text{in} \ (0,T] \ \times \ \Omega_{S} \\
 -\boldsymbol \nabla \cdot \left( \, \epsilon_{r}^{S} \, \boldsymbol \nabla \Phi \right)
\; &= \;  
\frac{q}{\epsilon_{0}} \left[ 
\rho_{p}^{e} - 
\rho_{n}^{e}  - (\rho_{n} -\rho_{p}) \right] 
  && \text{in} \ (0,T] \ \times \ \Omega_{S}  \nonumber .
 \end{align} \f]


Where \f$\alpha_{i} \f$, \f$\mu_{i}\f$ and \f$D_{i}\f$ are the charge numbers mobility and 
diffusivity of carrier \f$i = n,p\f$
respectively.  The functions \f$\rho_{n}^{e}\f$
and \f$\rho_{p}^{e}\f$ are the equilibrium electron and hole densities respectively.  
The constants
\f$q\f$ and \f$\epsilon_{0}\f$ charge of the electron and permittivity of free space. 
The material permittivity is \f$\epsilon_{r}\f$ is assumed to be an invertible tensor.


We use Shockley-Reed-Hall recombination SOLARCELL::SRH_Recombination as our sink functional,
\f[ \begin{equation}
R(\rho_{n}, \rho_{p})
\; = \;
\frac{\rho_{i}^{2} \, -\, \rho_{n} \, \rho_{p}}
{\tau_{n} \, (\rho_{n} \, + \, \rho_{i}) \, + \, \tau_{p} \, ( \rho_{p} \, + \, \rho_{i})}.
\end{equation}
\f]


The term \f$\rho_{i}\f$ is the intrinsic electron density and can be spatially varying. 
\f$\tau_{n}, \ \tau_{p} \f$ are constants called the electron and hole lifetimes.  
The generation of electrons and holes is modeled using a macroscopic source function, 
Generation,


\f[
\begin{equation}
G(\textbf{x}) \; = \; \left\{ 
\begin{array}{lr}
\alpha (\textbf{x}) 
\, G_{0} \, e^{- \, \int_{0}^{s} \, \alpha 
(\, \textbf{x}_{0} \, + \, s' \, \boldsymbol \theta_{0} \, ) \, ds'} 
\qquad & \qquad \text{if} \quad  \textbf{x} 
\; = \; \textbf{x}_{0} \, + \, s \, \boldsymbol \theta_{0} \\
0 \qquad & \qquad \text{otherwise}
\end{array} \right. 
\end{equation} \f]



The point \f$\textbf{x}_{0}\f$ is the photon's incident location and 
\f$\boldsymbol \theta_{0}\f$ is the incident direction. The absorption coefficient 
\f$\alpha(\textbf{x}) \f$ has been averaged over all 
energy values of light that generate free carriers.  The term 
\f$G(\textbf{x}_{0}) \; [ \, \text{cm}^{-2} \, \text{s}^{-1} \, ]\f$ 
represents the surface photon flux at the point \f$\textbf{x}_{0}\f$. 


The portion of the boundary of the semiconductor \f$\Gamma_{S,D}\f$ is 
an Ohmic metal contact where the charge densities take on their 
equilibrium values (Electrons_Equilibrium and Holes_Equilibrium),
\f[ \begin{align}
\rho_{n}   \; &= \;  
\rho_{n}^{e}  && \text{on} \ (0,T] \ \times \ \Gamma_{S,D}  \\
\rho_{p} \; &= \; 
\rho_{p}^{e}  && \text{on} \ (0,T] \ \times \ \Gamma_{S,D}
\end{align} 
\f]

The potential on an Ohmic contact is the sum of the applied voltage 
\f$\Phi_{\text{app}} \f$ Applied_Bias and the so-called ``built-in'' potential 
\f$\Phi_{\text{bi}}\f$ Built_In_Bias,

\f[
\begin{align} 
\Phi  \;  = \; \Phi_{\text{bi}} \, + \, \Phi_{\text{app.}}  
&& \text{on} \ (0,T] \ \times \ \Gamma_{S,D}.
\end{align}
\f]

The location of \f$ \Gamma_{S,D} \f$ is decided in Grid_Maker::Grid::make_Dirichlet_boundaries. 
The portion of the boundary of the semiconductor \f$\Gamma_{S,N}\f$ is 
an oxide where the normal component of currents and electric field are zero,

\f[ \begin{align}
\left(  - \alpha_{n} \, \mu_{n}  \,  \,\boldsymbol \nabla \Phi \, \rho_{n} 
\, - \,  D_{n}  \, 
\boldsymbol \nabla \, \rho_{n}  \, \right)  \ \cdot \ \boldsymbol \eta \; &= \; 0 
					&& \text{on} \ (0,T] \ \times \ \Gamma_{S,N} \\
\left(  - \alpha_{p} \, \mu_{p}  \,  \,\boldsymbol \nabla \Phi \, \rho_{p} 
\, - \,  D_{p}  \, 
\boldsymbol \nabla \, \rho_{n}  \, \right)  \ \cdot \ \boldsymbol \eta \; &= \; 0 
					&& \text{on} \ (0,T] \ \times \ \Gamma_{S,N} \\
-\boldsymbol \nabla \Phi  \ \cdot \ \boldsymbol \eta \; &= \; 0 
					&& \text{on} \ (0,T] \ \times \ \Gamma_{S,N} \\
\end{align}
\f]

The location of \f$ \Gamma_{S,N} \f$ is set in 
Grid_Maker::Grid::make_Neumann_boundaries.


In the electrolyte component \f$(\Omega_{E})\f$ the transport of reductants 
(\f$ \rho_{r} \f$) and oxidants (\f$ \rho_{o} \f$ ) is governed by a 
similar drift-diffusion-Poisson system,

\f[
\begin{align}
\frac{\partial \rho_{r}}{\partial t}
\, + \, 
\boldsymbol \nabla  \cdot \, \left( - \alpha_{r} \, \mu_{r}   \, \boldsymbol \nabla \Phi \, \rho_{r} 
\, - \, D_{r} \,
 \boldsymbol \nabla \, \rho_{r}  \, \right) 
\; &= \; 0 , && \text{in} \ (0,T] \ \times \ \Omega_{E}  \nonumber \\
\frac{\partial \rho_{o}}{\partial t} 
\, + \, 
\boldsymbol \nabla \cdot  \, \left( - \alpha_{o} \, \alpha_{o}  \, \boldsymbol \nabla \Phi \, \rho_{o} 
\, - \, D_{o} \,
\boldsymbol \nabla \, \rho_{o} \, \right) 
\; &= \; 0, && \text{in} \ (0,T] \ \times \ \Omega_{E}  \\
-\boldsymbol \nabla \cdot \left( \, \epsilon^{E}_{r} \, \boldsymbol \nabla \Phi \right)
\; &= \;
- \, \frac{q}{\epsilon_{0}} \, ( \rho_{r} - \rho_{o}). && \text{in} \ (0,T] \ \times \ \Omega_{E} \nonumber
\end{align}
\f]

Where \f$\alpha_{i} \f$, \f$\mu_{i}\f$ and \f$D_{i}\f$ are the charge numbers mobility
 and diffusivity of carrier \f$i = r,o\f$ respectively. The electrolyte permittivity is
\f$\epsilon_{r}^{E}\f$.  The lack of doping profile in our electrolyte reflects the fact 
that it is charge neutral.  Our model does not include any generation or recombination 
mechanisms in the electrolyte domain since we are only considering so-called 
``heterogenous reactions." That is chemical reactions can only occur at the interface 
and not within the bulk of the electrolyte.  

We assume the interface is isolated  and that the electrolyte variables take on their 
bulk values (Reductants_Equilibrium, Oxidants_Equilibrium, Bulk_Bias) 
on boundary \f$\Gamma_{E,D}\f$,

\f[
\begin{align}
\rho_{r} \vert_{\Gamma_{E,D}}  \; = \;  \rho_{r}^{\infty}   ,
&&
\rho_{o} \vert_{\Gamma_{E,D}} \; = \; \rho_{o}^{\infty}  ,
&&
\Phi \vert_{\Gamma_{E,D}} \; = \; \Phi^{\infty} , 
&& \text{on} \ (0,T] \ \times \ \Gamma_{E,D} 
\end{align}
\f]

The location of \f$\Gamma_{E,D}\f$ is again set in 
Grid_Maker::Grid::make_Dirichlet_boundaries.
The insulating portion of the electrolyte \f$\Gamma_{E,N} \f$ is set in 
Grid_Maker::Grid::make_Neumann_boundaries.  On the insulting portion of the 
electrolyte we have,

\f[ \begin{align}
\left(  - \alpha_{r} \, \mu_{r}  \,  \,\boldsymbol \nabla \Phi \, \rho_{r} 
\, - \,  D_{r}  \, 
\boldsymbol \nabla \, \rho_{n}  \, \right)  \ \cdot \ \boldsymbol \eta \; &= \; 0 
					&& \text{on} \ (0,T] \ \times \ \Gamma_{E,N} \\
\left(  -\alpha_{o} \, \mu_{o}  \,  \,\boldsymbol \nabla \Phi \, \rho_{o} 
\, - \,  D_{o}  \, 
\boldsymbol \nabla \, \rho_{n}  \, \right)  \ \cdot \ \boldsymbol \eta \; &= \; 0 
					&& \text{on} \ (0,T] \ \times \ \Gamma_{E,N} \\
-\boldsymbol \nabla \Phi  \ \cdot \ \boldsymbol \eta \; &= \; 0 
					&& \text{on} \ (0,T] \ \times \ \Gamma_{E,N} \\
\end{align}
\f]

The default values of the boundaries of the semiconductor and electrolyte domain are
 set to be the interface \f$\Sigma\f$, so any portion of the boundaries which are not set
 by the object Grid_Maker::Grid are defined as the interface. The potential and
 displacement electric field are continuous across \f$\Sigma\f$,

\f[
\begin{align}
\Phi \, \vert_{\Sigma^{-}} \, = \, \Phi \, \vert_{\Sigma^{+}} 
&&  -\epsilon_{r}^{S} \, \boldsymbol \nabla \, \Phi 
\ \cdot \ \boldsymbol \eta_{\Sigma} \ \vert_{\Sigma^{-}} 
\, = \,
-\epsilon^{E}_{r} \, \boldsymbol \nabla \, \Phi \,
\cdot \ \boldsymbol \eta_{\Sigma} \ \vert_{\Sigma^{+}} 
&& \text{on} \ (0,T] \ \times \ \Sigma
\end{align}
\f]

The chemical reactions of the charge carriers on the interface are modeled using 
the following boundary conditions on the currents,

\f[
\begin{align}
\textbf{J}_{n} \cdot \boldsymbol \eta_{\Sigma}
\; = \; 
\left( - \alpha_{n} \, \mu_{n}  \,  \,\boldsymbol \nabla \Phi \, \rho_{n} 
\, - \,  D_{n}  \, 
\boldsymbol \nabla \, \rho_{n}  \, \right)  \cdot \boldsymbol \eta_{\Sigma} 
\; &= \; 
k_{et} \, ( \, \rho_{n} \, - \, \rho_{n}^{e} \, ) \, \rho_{o} , \\ 
\textbf{J}_{p} \cdot \boldsymbol \eta_{\Sigma}
\; = \; 
\left( - \alpha_{p} \, \mu_{p}  \, \boldsymbol \nabla \Phi \, \rho_{p} 
\, - \, D_{p}  \,
\boldsymbol \nabla \, \rho_{p} \, \right) \cdot \boldsymbol \eta_{\Sigma} 
\; &= \; 
k_{ht} \, ( \, \rho_{p} \, - \, \rho_{p}^{e} \, ) \, \rho_{r},  \\ 
\textbf{J}_{r} \cdot \boldsymbol \eta_{\Sigma}
\; = \; 
\left(  - \alpha_{r} \, \mu_{r}  \,\boldsymbol \nabla \Phi \, \rho_{r} 
\, - \, 
D_{r}  \,
\boldsymbol \nabla \, \rho_{r}  \, \right) \cdot \boldsymbol \eta_{\Sigma} \;  
&=
\;  \, k_{ht} \,  (\rho_{p} \, - \, \rho_{p}^{e}) \rho_{r} \, - \, 
k_{et} \, (\rho_{n} \, - \, \rho_{n}^{e}) \, \rho_{o}  , \\
\textbf{J}_{o} \cdot \boldsymbol \eta_{\Sigma}
\; = \; 
\left( - \alpha_{o} \, \mu_{o} \  \boldsymbol \nabla \Phi \, \rho_{o} 
\, - \, 
D_{o} \ 
\boldsymbol \nabla \, \rho_{o} \, \right) \cdot \boldsymbol \eta_{\Sigma} 
\; &= \; 
- k_{ht} \, (\rho_{p} \, - \, \rho_{p}^{e}) \, \rho_{r} \, + \, 
k_{et} \, (\rho_{n} \, - \, \rho_{n}^{e}) \, \, \rho_{o}.
\end{align}
\f]


The initial conditions are taken to be,
\f[ \begin{align}
\rho_{n}   \; &= \;  
\rho_{n}^{e}  && \text{on} \ \{0\} \ \times \ \Omega_{S}  \\
\rho_{p} \; &= \; 
\rho_{p}^{e}  && \text{on} \ \{0\} \ \times \ \Omega_{S} \\
\rho_{r}   \; &= \;  \rho_{r}^{\infty}  && \text{on} \ \{0\} \ \times \ \Omega_{E}  \\
\rho_{o} \; &= \; \rho_{o}^{\infty}  && \text{on} \ \{0\} \ \times \ \Omega_{E}
\end{align} 
\f]

These functions are assigned as objects of the classes, 
Electrons_Equilibrium, Holes_Equilibrium, Reductants_Equilibrium, Oxidants_Equilibrium.


\note In scaling our equations we use Einstein's relations \f$ D \ = \ U_{T} \ \mu \f$ 
and singular perturbation scaling, for more information see the
<a href="http://www.ma.utexas.edu/users/gamba/papers/Semi-Elec-Interf.pdf">
paper</a> on this model.  This eliminates the need for inputing diffusion constants.

\note We take \f$ \alpha_{n} = -1\f$, \f$\alpha_{p} = 1\f$ and have the constraint, 
\f$ \alpha_{o} \, - \, \alpha_{r} = 1\f$.  
These can be set in the constructor, SOLARCELL::SolarCellProblem().

The output of these simulations will be the calculations of potential, electric field, 
charge densities and the current.  The although in our model we have eliminated the 
charge of an electron, \f$q\f$, the output current through the device will involve the
 charge of the electron and is defined to be:

\f[
\textbf{J}(\textbf{x}) \; = \;
\left\{
\begin{array}{cc}
-q \, \alpha_{n} \, \textbf{J}_{n} \, - \, q \, \alpha_{p} \, \textbf{J}_{p}, & 
\textbf{x} \in \Omega_{S} \\
-q \, \alpha_{r} \, \textbf{J}_{r} \, - \, q \, \alpha_{o} \, \textbf{J}_{o}, & 
\textbf{x} \in \Omega_{E} \\
\end{array}
\right.
\f]

\note This definition of the current is continuous across the interface so long as 
\f$ \alpha_{o} \, - \, \alpha_{r} = 1\f$.

\subsection Methods

The overall stragy is to create a domain decomposition where in each subdomain we have
a reaction-drift-diffusion system of equations for a pair of charge carriers.  The
potential and electric field live in a superset of these two domains and information is
passed between the two subdomains as well as the superset domain.


The drift-diffusion transport equations are numerically approximated by 
the a local discontinuous Galerkin
method in space, see LDG_System and the classes 
within this namesapace for more details.
Poisson's equation for the potential and electric field is approximated using a 
mixed finite element method see MixedPoisson::MixedFEM() for more details.


Time stepping is handeled in a specific way such that nonlinear terms are
linearized by time lagging.  Solutions to Poisson equations are updated using
implicit density values, while the charge carriers use an IMEX strategy. 
The overall time stepping strategy
is termed a "parallel Gummel-Schwarz method."  The steps in this algorithm are
presented in the flow chart, \image html parallel_gummel_schwarz.jpg
For more deail see LDG_System::LDG. 