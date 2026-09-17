#include "PbPbH011InjectionToys.C"

void PbPbH012InjectionToys(
    const char* key, const char* templateType, const char* templateEventWeight,
    const char* dataPath, const char* dataTree, const char* referencePath,
    const char* referenceTree, const char* selection,
    double yieldMinus, double yieldCentral, double yieldPlus,
    int toysPerEnsemble, int seedBase, const char* outputDirectory,
    bool runAsimov = true, bool runToys = true)
{
    PbPbH011InjectionToys(
        key, templateType, dataPath, dataTree, referencePath, referenceTree,
        selection, yieldMinus, yieldCentral, yieldPlus, toysPerEnsemble,
        seedBase, outputDirectory, runAsimov, runToys,
        templateEventWeight, 3.94, 28);
}
