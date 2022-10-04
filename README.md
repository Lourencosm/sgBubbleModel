# sgBubbleModel

sgBubbleModel is a sub-grid bubble model library implemented in OpenFOAM to be coupled with the interFoam solver.

# interFoam_with_sgbm

interFoam solver with a few addional entries to load and run the sgBubbleModel.

# turbModels

Variations of the k-epsilon and k-omega SST turbulence models that include bubble source/sink terms.

\*All libraries are valid for OpenFOAM-v1906 version.
An update for OpenFOAM-v2112 will be upload in late 2022.

## Author

Lourenço Sassetti da Silva Mendes, PhD.

[LNEC](http://www.lnec.pt/hidraulica-ambiente/en/team/lourenco-sassetti-mendes/)

[GIT](https://github.com/Lourencosm)

[RESEARCH GATE](https://www.researchgate.net/profile/Lourenco-Sassetti-Mendes)


This work was developed under the guidance of:

Javier L. Lara (https://ihcantabria.com/en/directorio-personal/responsable-de-grupo/javier-lopez-lara/)

Teresa Viseu (http://www.lnec.pt/hidraulica-ambiente/pt/equipa/maria-teresa-viseu/)


## Publications

Lourenço Sassetti Mendes, Javier L. Lara, and Maria Teresa Viseu. ‘Is the Volume-of-Fluid
Method Coupled with a Sub-Grid Bubble Equation Efficient for Simulating Local and
Continuum Aeration?’ In: Water 13.11 (2021), p. 1535. doi: 10.3390/w13111535. https://www.mdpi.com/2073-4441/13/11/1535

Lourenço Sassetti Mendes. ‘Computational Fluid Dynamics of Aerated Flows in Spillways
Chutes’. PhD. Santander, España: Universidad de Cantabria, 2021. http://hdl.handle.net/10902/23909


## Installation

1.  Load OpenFOAM bashrc. Depends on your OpenFOAM installation setup. Example:

```
source ~/OpenFOAM/OpenFOAM-v1906/etc/bashrc
```

2.  Compile the **sgBubbleModel**:

```
cd sgBubbleModel
./make_sgbm
```

3.  Compile the boundary condition **variableHeightFlowRateInletVelocityDoublePowerLawBL**:

```
cd boundaryConditions/variableHeightFlowRateInletVelocityDoublePowerLawBL_00/
./make
```

2.  Compile the **interFoam_with_sgbm**:

```
cd interFoam_with_sgbm
./make_interFoam_with_sgbm
```

## Usage

Run the **tutorials cases**:

```
cd tutorials/impingingJet
interFoam_with_sgbm
```
and
```
cd  tutorials/spillwayChute
interFoam_with_sgbm
```


## Contributing

1.  Fork it!
2.  Create your feature branch: `git checkout -b my-new-feature`
3.  Commit your changes: `git commit -am 'Add some feature'`
4.  Push to the branch: `git push origin my-new-feature`
5.  Submit a pull request :D

## Credits

[FCT](https://fct.pt/)
[IH Cantabria](https://ihcantabria.com)
[LNEC](http://www.lnec.pt/)

## License

Licensed under the GNU General Public License v3.0 - see the LICENSE.md file for details