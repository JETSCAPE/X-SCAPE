# Contributing

Contributions to JETSCAPE (new module, feature, bug fix, etc.) are welcome. To contribute, read through the instructions below and open a [Pull Request](https://github.com/JETSCAPE/JETSCAPE/pulls) with your changes, or an [Issue](https://github.com/JETSCAPE/JETSCAPE/issues) describing what you intend to do.

## Developing modules

To develop a new JETSCAPE module, you should inherit from the relevant base class (InitialState, JetEnergyLoss, etc.)
and implement the relevant initialization and execution functions, described in detail in [The JETSCAPE framework](https://arxiv.org/abs/1903.07706)
Section 3.3.

Additionally, you must register your module with the framework with the following steps:
- Add the following to your module .h:
  ```
  private:
  // Allows the registration of the module so that it is available to be used by the Jetscape framework.
  static RegisterJetScapeModule<MyClass> reg;
  ```
- Add the following to your module .cc:
  ```
  // Register the module with the base class
  RegisterJetScapeModule<MyClass> MyClass::reg("CustomModuleBlahBlah");
  ```
where `MyClass` is the name of your class, and "CustomModuleBlahBlah" is the name that should be added to the XML configuration.
You can see any of the established modules, e.g.  `Matter`, as an example.

Important Note: In the case of custom modules, you *must* start your module name with "CustomModule..."
in order for it to be recognized by the framework (for custom writers, you must start the name with "CustomWriter").

New modules should not use multiple inheritance, if avoidable.

Once these steps are done, one can just add the module name to the XML, and it will be automatically available to run in JETSCAPE.

## Git Management

Tips for git management are found on the corresponding [wiki page](https://github.com/JETSCAPE/JETSCAPE/wiki/Tips-for-git-management).

## Doxygen documentation

An overview of the code structure and functionality can be obtained by building the HTML doxygen documentation. Run
```
doxygen JetScapeDoxy.conf
```
in the JETSCAPE directory. Open the `index.html` file in the newly-created `html` directory in your favorite browser to see the documentation.
