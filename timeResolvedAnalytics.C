#include "timeResolvedAnalytics.H"

namespace Foam
{

timeResolvedAnalytics::timeResolvedAnalytics(const fvMesh& mesh, const dictionary& dict)
:
    mesh_(mesh),
    dict_(dict),
    timeHistory_(),
    temperatureHistory_(),
    pressureHistory_()
{}

void timeResolvedAnalytics::calculate(const volScalarField& T, const volScalarField& p)
{
    scalar currentTime = mesh_.time().value();
    scalar avgTemp = gAverage(T);
    scalar avgPressure = gAverage(p);

    timeHistory_.append(currentTime);
    temperatureHistory_.append(avgTemp);
    pressureHistory_.append(avgPressure);
}

void timeResolvedAnalytics::write()
{
    // Write time series data to file
    OFstream timeSeriesFile(mesh_.time().path()/"timeSeriesData.csv");
    
    timeSeriesFile << "Time,Temperature,Pressure" << nl;
    forAll(timeHistory_, i)
    {
        timeSeriesFile << timeHistory_[i] << ","
                       << temperatureHistory_[i] << ","
                       << pressureHistory_[i] << nl;
    }
}

}
