---
title: Contributed
permalink: /Contributed/
---

Some people extended SUMO or built tools to make it more usable. Not all
of these extensions have found its way to what we would call "SUMO
core".

# Within SUMO

The following extensions became a core part of the SUMO package.

- **[TraCI](TraCI.md)**
  
    online interaction with the simulation

# Included in the Distribution

The following contributions are included in the package, but are less
supported.

- **[Contributed/Cadyts](Contributed/Cadyts.md)**

    a tool by Gunnar Flötteröd which adapts the simulate flows to the
    real flows in a known net

- **[SUMOPy](Contributed/SUMOPy.md)**

    a tool by Joerg Schweizer supporting the whole SUMO toolchain with a
    GUI especially for demand modelling

- **[LiSuM](Tools/LiSuM.md)**

    a middleware that couples LISA+ and SUMO to simulate real-world
    traffic light controllers.

- ~~[Contributed/SUMO Traffic
Modeler](Contributed/SUMO_Traffic_Modeler.md)~~

    ~~a graphical user interface for defining high-level traffic for
    SUMO~~ (obsolete)

# External

The following extensions are managed and supported by other parties.

## Demand Generators

- **[Citymob for Roadmaps](http://www.grc.upv.es/Software/c4r.html)**

    mobility pattern generator for vehicular networks (VANet focus)

## Scenario and Network Editors

- **[eWorld](http://eworld.sourceforge.net/)**

    an application that allows to convert and enrich roads networks

- **[Sumo2Unreal](https://github.com/AugmentedDesignLab/Sumo2Unreal)**

    an importer for SUMO's .net.xml road network files into the Unreal
    Engine 4 environment.

## Connections to Network Simulators

- **[Veins](http://www7.informatik.uni-erlangen.de/veins/)**

    connects SUMO to OMNeT++

- **[<font color="#0174DF">Tra</font><font color="#FF0000">NS</font>](http://trans.epfl.ch/)**

    connects SUMO to ns-2

- **[MOVE](http://lens1.csie.ncku.edu.tw/wiki/doku.php?id=%E2%80%A7realistic_mobility_generator_for_vehicular_networks)**

    connects SUMO to ns-2

- **[VSimRTI](http://www.dcaiti.tu-berlin.de/research/simulation/)**

    connects SUMO to OMNeT++ and JiST/SWANS

## Other

- **[FLOW](https://flow-project.github.io/)**

    a framework for applying reinforcement learning and custom
    controllers to SUMO, developed

# Purgatory

The following extensions exist or have existed, but their state is
unclear.

- **[Contributed/iTranSIM](Contributed/iTranSIM.md)**

    extension by online-calibration using induction loop data by Tino
    Morenz

- **[Contributed/SmallMaps](Contributed/SmallMaps.md)**

    prunes road networks to a given boundary; by Andreas Florides

- **[Contributed/SUMOPlayer](Contributed/SUMOPlayer.md)**

!!! note
    SUMOPlayer was removed in release 0.24.0. You should be able to use [traceExporter.py](Tools/TraceExporter.md) for the same task.