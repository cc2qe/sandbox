

// returns true if testValue matches the flagToCheckFor 
bool testFlag(int testValue, int flagToCheckFor) {
    return (testValue & flagToCheckFor) == flagToCheckFor;
}
