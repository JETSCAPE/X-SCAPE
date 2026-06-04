# References

A consolidated bibliography for X-SCAPE and the physics models it integrates.
Please cite the relevant papers when you use a module for scientific work, and
**always cite the JETSCAPE framework paper**.

## Framework

| Topic | Reference |
|---|---|
| **JETSCAPE framework** (cite this) | The JETSCAPE Collaboration, *The JETSCAPE framework*, [arXiv:1903.07706](https://arxiv.org/abs/1903.07706) |
| X-SCAPE small-system multistage framework | *A multistage framework for the evolution of jets and high-$p_T$ probes in small collision systems*, [arXiv:2308.02650](https://arxiv.org/abs/2308.02650) |
| Soft–hard framework, exact four-momentum conservation | *A soft–hard framework with exact four-momentum conservation for small systems*, [arXiv:2407.17443](https://arxiv.org/abs/2407.17443) |
| Multistage Monte-Carlo jet modification | Cao & Majumder, [arXiv:1705.00050](https://arxiv.org/abs/1705.00050) |
| Multistage jet quenching (results) | JETSCAPE, [arXiv:2002.07124](https://arxiv.org/abs/2002.07124) |
| Constraints on jet quenching (Bayesian) | JETSCAPE, [arXiv:2009.02410](https://arxiv.org/abs/2009.02410) |

## Initial state

| Model | Reference |
|---|---|
| **TRENTo** | Moreland, Bernhard, Bass, [arXiv:1412.4708](https://arxiv.org/abs/1412.4708) |
| **3DGlauber** (dynamical 3-D initial state) | Shen & Schenke, [arXiv:1710.00881](https://arxiv.org/abs/1710.00881); [arXiv:1711.10544](https://arxiv.org/abs/1711.10544) |
| **IP-Glasma** | Schenke, Tribedy, Venugopalan, [arXiv:1202.6646](https://arxiv.org/abs/1202.6646) |
| **SMASH initial condition** | Weil et al., [arXiv:1606.06642](https://arxiv.org/abs/1606.06642) |

## Pre-equilibrium

| Model | Reference |
|---|---|
| Free-streaming + Landau matching (freestream-milne) | free-streaming pre-equilibrium approach; see the `freestream-milne` package documentation |
| Glasma classical Yang–Mills | Schenke, Tribedy, Venugopalan, [arXiv:1202.6646](https://arxiv.org/abs/1202.6646) |

## Hydrodynamics

| Model | Reference |
|---|---|
| **MUSIC** | Schenke, Jeon, Gale, [arXiv:1004.1408](https://arxiv.org/abs/1004.1408); Paquet et al., [arXiv:1509.06738](https://arxiv.org/abs/1509.06738) |
| **CLVisc** | Pang, Petersen, Wang, [arXiv:1802.04449](https://arxiv.org/abs/1802.04449) |
| Gubser flow (analytic benchmark) | Gubser, [arXiv:1006.0006](https://arxiv.org/abs/1006.0006) |

## Hard process

| Model | Reference |
|---|---|
| **PYTHIA 8** | Sjöstrand et al., [arXiv:1410.3012](https://arxiv.org/abs/1410.3012) |

## Jet energy loss

| Model | Reference |
|---|---|
| **MATTER** | Majumder, [arXiv:1301.5323](https://arxiv.org/abs/1301.5323); Cao & Majumder, [arXiv:1712.10055](https://arxiv.org/abs/1712.10055) |
| **LBT** (Linear Boltzmann Transport) | He, Luo, Wang, Zhu, [arXiv:1503.03313](https://arxiv.org/abs/1503.03313) |
| **MARTINI** | Schenke, Gale, Jeon, [arXiv:0909.2037](https://arxiv.org/abs/0909.2037) |
| **AdS/CFT** holographic energy loss | Chesler & Rajagopal, [arXiv:1402.6756](https://arxiv.org/abs/1402.6756); Casalderrey-Solana et al., [arXiv:1101.0618](https://arxiv.org/abs/1101.0618) |

## Hadronization

| Model | Reference |
|---|---|
| Lund string (colorless / colored) | PYTHIA 8 — Sjöstrand et al., [arXiv:1410.3012](https://arxiv.org/abs/1410.3012) |
| **Hybrid hadronization** (recombination + string) | Han, Fries, Cao, [arXiv:1601.01583](https://arxiv.org/abs/1601.01583) |

## Particlization & afterburner

| Model | Reference |
|---|---|
| **iSS** Cooper–Frye sampler | Shen, Qiu, Song, Bernhard, Bass, Heinz (iEBE-VISHNU), [arXiv:1409.8164](https://arxiv.org/abs/1409.8164) |
| **Cornelius** freeze-out surface finder | Huovinen & Petersen, [arXiv:1206.3371](https://arxiv.org/abs/1206.3371) |
| **SMASH** hadronic transport | Weil et al., [arXiv:1606.06642](https://arxiv.org/abs/1606.06642); [arXiv:1808.06832](https://arxiv.org/abs/1808.06832) |

## Output formats

| Topic | Reference |
|---|---|
| HepMC heavy-ion conventions | [arXiv:1912.08005](https://arxiv.org/abs/1912.08005) |

## Tunes and further material

- **Published parameter tunes:** [JETSCAPE/Default-tunes](https://github.com/JETSCAPE/Default-tunes)
- **Summer-school material:** repositories under the [JETSCAPE organization](https://github.com/JETSCAPE)
- **Project site:** [jetscape.org](http://jetscape.org)

!!! note "Accuracy of references"
    The arXiv identifiers above are the canonical community references for each
    model. For a few models (e.g. the specific AdS/CFT prescription and the
    freestream-milne implementation as wired into X-SCAPE) consult the source
    headers in `src/` and the upstream package repositories for the precise
    citation matching the version you build, and prefer the citation requested
    by the upstream package where one is given.
