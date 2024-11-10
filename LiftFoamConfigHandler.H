#ifndef LiftFoamConfigHandler_H
#define LiftFoamConfigHandler_H

#include "dictionary.H"
#include "IOdictionary.H"
#include "IOstreams.H"
#include "fvCFD.H"
#include "volFields.H"

namespace Foam
{

class LiftFoamConfigHandler
{
public:
    // Constructor
    LiftFoamConfigHandler(const Time& runTime, const fileName& dictPath);

    // Access functions
    const dictionary& dict() const { return dict_; }
    const Time& runTime() const { return runTime_; }

    // Configuration validation
    bool validateConfiguration() const;
    bool validateFields() const;
    
    // Entry retrieval with validation
    template<class T>
    T getEntry(const word& entryName, const T& defaultValue) const
    {
        if (dict_.found(entryName))
        {
            return dict_.get<T>(entryName);
        }
        else if (dict_.found("LIFT") && dict_.subDict("LIFT").found(entryName))
        {
            return dict_.subDict("LIFT").get<T>(entryName);
        }
        
        Info<< "Using default value for " << entryName << ": " << defaultValue << endl;
        return defaultValue;
    }

    // Required entry retrieval
    template<class T>
    T getRequiredEntry(const word& entryName) const
    {
        if (dict_.found(entryName))
        {
            return dict_.get<T>(entryName);
        }
        else if (dict_.found("LIFT") && dict_.subDict("LIFT").found(entryName))
        {
            return dict_.subDict("LIFT").get<T>(entryName);
        }
        
        FatalErrorInFunction
            << "Required entry '" << entryName << "' not found in configuration."
            << abort(FatalError);
            
        return T();
    }

private:
    const Time& runTime_;
    const IOdictionary dict_;

    // Validation helpers
    template<class T>
    bool checkEntry(const word& entryName, const T& defaultValue) const
    {
        const T value = getEntry<T>(entryName, defaultValue);
        if (value == defaultValue)
        {
            Info<< "Warning: Entry '" << entryName 
                << "' using default value: " << defaultValue << endl;
        }
        return true;
    }

    bool validateLaserParameters() const;
    bool validateMaterialParameters() const;
    bool validateNumericalParameters() const;
    bool validateFieldDimensions(const fvMesh& mesh) const;
};

}

#endif
