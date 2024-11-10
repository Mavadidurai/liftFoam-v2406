#ifndef ValidationReporter_H
#define ValidationReporter_H

#include "fvCFD.H"
#include "dictionary.H"
#include <map>
#include <string>

namespace Foam
{

class ValidationReporter
{
private:
    struct ValidationResult
    {
        bool passed;
        scalar value;
        scalar limit;
        std::string message;
        word field;
        
        ValidationResult()
        :
            passed(true),
            value(0),
            limit(0),
            message(""),
            field("")
        {}
    };

    const Time& runTime_;
    const dictionary& dict_;
    std::map<word, ValidationResult> results_;
    label writeFrequency_;
    bool stopOnFailure_;
    bool warnOnly_;
    label logLevel_;
    bool writeFields_;
    
    void writeValidationReport() const;
    void writeValidationFields() const;

public:
    ValidationReporter(const Time& runTime, const dictionary& dict);

    void addResult
    (
        const word& name,
        bool passed,
        scalar value,
        scalar limit,
        const std::string& message,
        const word& field
    );

    void clearResults();
    bool allPassed() const;
    void report();
    
    static word validationResultsFile;
};

} // End namespace Foam

#endif

