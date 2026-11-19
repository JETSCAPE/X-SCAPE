# Colour flow through the jet energy loss, and iColored hadronization with LBT

This document describes how the partons of Matter and LBT carry colour tags, and how
`iColoredHadronization` closes the colour lines that scatterings off the medium leave
open, so that a Matter + LBT shower can be hadronized through its colour flow (Pythia
strings), as the O+O ISR chain does. Branch `lbt_colored_hadro`; design, checks and
validation in js-contrib `docs/Plans/PLAN_lbt_colour_tags.md`.

## What changed

* **One colour-tag counter per event.** `HardProcess::NextColorTag()` hands out tags from
  a counter that is reset every event and starts above every tag the hard process and
  iMATTER use (Pythia's own tags, the per-parton blocks of 1000, iMATTER's ISR blocks).
  Matter and LBT draw their new tags from it. Before, Matter's counter was a member of the
  Matter object: the energy-loss task of the first shower lives for the whole run, its
  counter only grows, and the other showers' Matter objects are copies of it made each
  event. Once that counter was above the per-parton blocks, every shower numbered its tags
  from the same value. It got there through the ISR siblings' bases (which grow over the
  job) and through the uninitialized `max_color` Matter wrote into the daughters of a photon
  emission. In O+O 131 of 150 events had a tag used more than once, and so did 83 % of
  the events of a pp Z+jet run (PythiaGun + Matter); Pythia strings them without a
  message. Details in the plan, "Where the duplicates came from".
* **Matter's recoils** (`recoil_on` 1) keep the large-Nc colour flow `Matter::flavor`
  computes: the recoil its outgoing tags, the continuing jet parton its outgoing tags, the
  hole (pstat −1) the tags of the incoming thermal parton. Before, recoil and hole were set
  colourless.
* **LBT** keeps the tags of the partons it receives (before, every parton lost them at its
  first LBT step), gives every elastic scattering the flow of its channel (Matter's rules;
  channels 2, 6 and 8 added) and every radiated gluon Matter's splitting rule.
* **iColoredHadronization**
  * hands holes (pstat −1) and the liquefier's partons (−11, −13, −17) to Pythia no more;
    the holes' energy is recorded. The beam remnants always go to Pythia;
  * fails an event with duplicate tags (or relabels them, `<colour_duplicates>`);
  * closes the open colour lines (`<colour_closure>`, below);
  * checks Pythia: a fragmentation that fails is retried once with soft end partons and
    otherwise throws (before, the input partons were written out as "hadrons");
  * logs a colour summary at the end of the run.

No random number is drawn for any tag: the kinematics of every chain are unchanged, and
ColorlessHadronization (which does not read the tags) gives the same hadrons as before.
The hadrons of iColoredHadronization do change (no duplicates any more): runs before and
after are not combinable.

## The colour lines that end in the medium

A scattering off a thermal parton exchanges a colour line with the medium. The thermal
parton (the hole) is not a string parton, so every such scattering leaves one line whose end
is in the medium: an open colour (a quark hole), an open anticolour (an antiquark hole) or
one of each (a gluon hole). `<JetHadronization><colour_closure>` decides how they are
closed:

| `<colour_closure>` | what it does |
|---|---|
| `reconnect` (default) | join open colour and anticolour ends within a shower, the shortest string (smallest invariant mass) first; insert the holes' medium partners (the conjugate parton with the hole's tags swapped and its momentum) for what is left |
| `reconnect_global` | the same, joining across showers too (a systematic: string drag between the jets) |
| `thermal` | every open line gets its hole's medium partner (the literal colour flow) |
| `fake` | a soft (1 GeV) light (anti)quark along the beam for every open end |
| `none` | throw on an open line (the behaviour before) |

With every mode but `fake` and `none`, a line still open after the partners (broken without
a hole of its own, e.g. by the liquefier) is joined across showers, and what is then left
gets a soft end parton; both are counted.

## XML

```xml
<JetHadronization>
  <name>iColored</name>
  <colour_closure>reconnect</colour_closure>    <!-- reconnect | reconnect_global | thermal | fake | none -->
  <colour_duplicates>fail</colour_duplicates>   <!-- fail (the event throws) | relabel (a counted repair) -->
  <reseed_per_event>0</reseed_per_event>        <!-- 1: reseed Pythia every event (reproducible per event) -->
  <store_input>0</store_input>                  <!-- 1: keep the parton list handed to Pythia (PyJetscape) -->
</JetHadronization>
```

Defaults are in `config/jetscape_main.xml`. Matter (`recoil_on`) and LBT need no switch for
their tags.

### Matter alone: only the tag fix

With Matter alone (`recoil_on` 0, no `<Lbt>`) there are no holes, so no colour line ends in
the medium and `<colour_closure>` never acts. Against the code before these changes, only the
one counter per event (no duplicate tags) and iColored's checks (a failed fragmentation is no
longer written out as "hadrons") change the jet hadrons. This is the reference for studies of
the closure or of LBT's tags (Matter with recoils, Matter + LBT): compare against it, not
against runs made before the fix. `none` makes it explicit, since an open line then throws
instead of being closed:

```xml
<Eloss>
  <Matter> ... <recoil_on>0</recoil_on> ... </Matter>   <!-- and no <Lbt> -->
</Eloss>
<JetHadronization>
  <name>iColored</name>
  <colour_closure>none</colour_closure>
</JetHadronization>
```

The colour summary at the end of the run then shows no holes, joins, partners or soft end
partons. In js-contrib's O+O job this is `--matter-only --colour-closure none`.

## Running it with runJetscape

`config/jetscape_user_OO_ISR_MatterLBT_iColored.xml`: O+O 200 GeV, ISR (PythiaIsrGun +
iMATTER), MUSIC, Matter (`recoil_on` 1) + LBT, iColoredHadronization, final partons and
hadrons written as ASCII. From the X-SCAPE source directory, in a working directory of its
own (see `README_Launcher.md`):

```bash
./examples/run_in_workdir.sh -d build -w /tmp/oo_lbt \
    -u config/jetscape_user_OO_ISR_MatterLBT_iColored.xml -m config/jetscape_main.xml
```

The log ends with the run's colour summary, e.g.

```
iColoredHadronization: colour summary (closure reconnect): 5 events hadronized, 0 failed (thrown); 5 with open colour lines
iColoredHadronization:   holes 33 (25.2 GeV, not hadronized), open ends 44, joined pairs 13 (across showers as a fallback 0), medium partners 14 (11.1 GeV), soft end partons 0 (0 GeV), retries 0, duplicate tags 0, liquefier partons skipped 0
```

"Failed" events, joins across showers, soft end partons and retries should be rare; if they
are not, the colour flow is broken somewhere.

**Things to know**

* **A failing event stops runJetscape.** iColoredHadronization throws on duplicate tags
  (with `fail`), on an open line with `none`, and when Pythia fails also with soft end
  partons, and runJetscape catches no exception. PyJetscape's O+O driver skips such an event
  and records it.
* **The per-event audit** (what each event contained and how it was closed) is available
  from C++ (`iColoredHadronization::GetLastColourAudit`, `GetColourTotals`) and from
  PyJetscape; runJetscape writes only the summary above.
* **The holes are not hadronized**: the jet hadrons are higher by the holes' energy minus
  the partners' (both in the summary). Negative hadrons for the holes, as
  HybridHadronization makes them, are not implemented.
* **With the liquefier**, the partons it absorbs take their tags with them; those lines are
  closed by the fallbacks. Not validated yet.

## From Python (PyJetscape)

js-contrib `contribs/PyJetscape/README.md`, "Jet hadronization with colour
(iColoredHadronization)": the per-event audit, choosing the closure, hadronizing stored
partons with iColored. js-contrib's O+O production (`contribs/PyJetscape/example/
prod_OO_isr/`, `config/OO_MCGlauber_MUSIC_isr.xml`) runs Matter + LBT with iColored by
default (`--matter-only`: Matter alone, as before).

## Related work

| here | closest counterpart elsewhere | reference |
|---|---|---|
| colour of the scatterings off the medium (the holes' open lines) | colour-differential energy loss: colour flows from the leading parton to the medium in most medium-modified branchings, giving heavier singlets and softer hadrons | Beraudo, Milhano, Wiedemann, Phys. Rev. C 85 (2012) 031901, arXiv:1109.5025; JHEP 07 (2012) 144, arXiv:1204.4342 |
| `reconnect` | JEWEL: recoils turned into gluons and inserted into the jet's string next to the parton they scattered with | Kunnawalkam Elayavalli, Zapp, JHEP 07 (2017) 141, arXiv:1707.01539 |
| joining by the shortest string | colour reconnection in Pythia (string length λ) and Herwig (cluster masses) | Andersson, Gustafson, Söderberg, Nucl. Phys. B 264 (1986) 29; Sjöstrand, Khoze, Z. Phys. C 62 (1994) 281, hep-ph/9310242; Christiansen, Skands, JHEP 08 (2015) 003, arXiv:1505.01681; Gieseke, Röhr, Siódmok, Eur. Phys. J. C 72 (2012) 2225, arXiv:1206.0041 |
| `reconnect_global` | Pythia/Angantyr global colour reconnection in heavy-ion events | Bierlich, Gustafson, Lönnblad, Shah, JHEP 10 (2018) 134, arXiv:1806.10820; Eur. Phys. J. C 84 (2024) 231, arXiv:2309.12452 |
| `thermal` | Hybrid Hadronization: open lines completed with thermal partons from the hydro hypersurface, with recombination; POWLANG: heavy quarks string to thermal light partons | Han, Fries, Ko, Phys. Rev. C 93 (2016) 045207, arXiv:1601.00708; JETSCAPE (Sengupta et al.), arXiv:2501.16482; Beraudo et al., Eur. Phys. J. C 82 (2022) 607, arXiv:2202.08732 |
| `fake` | ColoredHadronization's fake beam partons | JETSCAPE (Kumar et al.), Phys. Rev. C 102 (2020) 054906, arXiv:1910.05481; Putschke et al., arXiv:1903.07706 |
| holes not hadronized | LBT's negative partons; JEWEL's thermal-momentum subtraction | He, Luo, Wang, Zhu, Phys. Rev. C 91 (2015) 054908, arXiv:1503.03313; Milhano, Wiedemann, Zapp, Phys. Lett. B 779 (2018) 409, arXiv:1707.04142 |

The ISR chain these hadronize is that of Soudi, Zhao, Majumder et al. (JETSCAPE), Phys. Rev. C
112 (2025) 014905, arXiv:2407.17443. Differences to the counterparts and the full reference
list: js-contrib `docs/Plans/PLAN_lbt_colour_tags.md`, "Related work and references".

## Validation (GB10, 2026-10-21)

Two 75-event O+O jobs per setting: Matter only (identical partons, no duplicate tag any
more), Matter with recoils (77 holes) and Matter + LBT (632 holes): every line closed by
joins and partners only (no soft end parton, no cross-shower join, no retry, no failure),
every tag paired exactly once in what went to Pythia, no parton among the hadrons, hadron
energy = Pythia input energy; `thermal` closes everything with the holes' partners alone,
which checks the flow of every channel. The Au+Au two-stage production (Matter with
recoils + LBT + liquefier, ColorlessHadronization) is identical except the stored colour
tags.
