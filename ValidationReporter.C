#include "ValidationReporter.H"
#include "fvMesh.H"
#include "volFields.H"
#include "OFstream.H"
#include <string>
#include <cstdlib>

namespace Foam {

word ValidationReporter::validationResultsFile = "validationResults";

ValidationReporter::ValidationReporter
(
    const Time& runTime,
    const dictionary& dict
)
:
    runTime_(runTime),
    dict_(dict),
    writeFrequency_
    (
        dict.subDict("validationControls").get<label>("writeFrequency")
    ),
    stopOnFailure_
    (
        dict.subDict("validationControls").get<bool>("stopOnFailure")
    ),
    warnOnly_
    (
        dict.subDict("validationControls").get<bool>("warnOnly")
    ),
    logLevel_
    (
        dict.subDict("validationControls").get<label>("logLevel")
    ),
    writeFields_
    (
        dict.subDict("validationControls").get<bool>("writeFields")
    )
{}

void ValidationReporter::addResult
(
    const word& name,
    bool passed,
    scalar value,
    scalar limit,
    const std::string& message,
    const word& field
)
{
    ValidationResult result;
    result.passed = passed;
    result.value = value;
    result.limit = limit;
    result.message = message;
    result.field = field;
    
    results_[name] = result;

    if (!passed && stopOnFailure_ && !warnOnly_)
    {
        FatalErrorInFunction
            << "Validation failed: " << name << nl
            << message << nl
            << "Value: " << value << nl
            << "Limit: " << limit << nl
            << "Field: " << field
            << abort(FatalError);
    }
}

void ValidationReporter::clearResults()
{
    results_.clear();
}

bool ValidationReporter::allPassed() const
{
    for (const auto& result : results_)
    {
        if (!result.second.passed)
        {
            return false;
        }
    }
    return true;
}

void ValidationReporter::report()
{
    if (runTime_.timeIndex() % writeFrequency_ == 0)
    {
        writeValidationReport();
        if (writeFields_)
        {
            writeValidationFields();
        }
    }
}

void ValidationReporter::writeValidationReport() const
{
    fileName reportFile
    (
        runTime_.path()/"postProcessing"/validationResultsFile/
        runTime_.timeName()/"validationReport.txt"
    );

    mkDir(reportFile.path());
    
    OFstream os(reportFile);
    
    os  << "Validation Report - Time: " << runTime_.timeName() << nl
        << "----------------------------------------" << nl;

    for (const auto& result : results_)
    {
        os  << "Test: " << result.first << nl
            << "  Passed: " << (result.second.passed ? "YES" : "NO") << nl
            << "  Value: " << result.second.value << nl
            << "  Limit: " << result.second.limit << nl
            << "  Message: " << result.second.message << nl
            << "  Field: " << result.second.field << nl
            << "----------------------------------------" << nl;
    }
}

void ValidationReporter::writeValidationFields() const
{
    fileName fieldsDir
    (
        runTime_.path()/"postProcessing"/validationResultsFile/
        runTime_.timeName()/"fields"
    );
    mkDir(fieldsDir);

    for (const auto& result : results_)
    {
        if (!result.second.passed && result.second.field != "")
        {
            const fvMesh& mesh = 
                runTime_.lookupObject<fvMesh>(polyMesh::defaultRegion);
                
            if (mesh.foundObject<volScalarField>(result.second.field))
            {
                const volScalarField& field = 
                    mesh.lookupObject<volScalarField>(result.second.field);
                    
                field.write();
            }
        }
    }
}

} // End namespace Foam

