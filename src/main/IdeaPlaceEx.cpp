#include "IdeaPlaceEx.h"
/* Parsers */
#include "parser/ProgArgs.h"
#include "parser/ParserTechSimple.h"
#include "parser/ParserPin.h"
#include "parser/ParserConnection.h"
#include "parser/ParserNetwgt.h"
#include "parser/ParserGds.h"
#include "parser/ParserSymFile.h"

PROJECT_NAMESPACE_BEGIN

// A global pointer to NLP
// This is for a trick to overly use static members in the NlpWnconj class
NlpWnconj *nlpPtr = nullptr;

void IdeaPlaceEx::readTechSimpleFile(const std::string &techsimple)
{
    ParserTechSimple(_db).read(techsimple);
}

void IdeaPlaceEx::readPinFile(const std::string &pinFile)
{
    ParserPin(_db).read(pinFile);
}

void IdeaPlaceEx::readConnectionFile(const std::string &connectionFile)
{
    ParserConnection(_db).read(connectionFile);

    // Init cells before read in the gds
    if(!_db.initCells())
    {
        ERR("IdeaPlaceEx::%s initializing the cells failed! \n", __FUNCTION__);
        Assert(false);
    }
}

void IdeaPlaceEx::readNetWgtFile(const std::string &netWgtFile)
{
    ParserNetwgt(_db).read(netWgtFile);
}

void IdeaPlaceEx::readSymFile(const std::string &symFile)
{
    ParserSymFile(_db).read(symFile);
}
void IdeaPlaceEx::readGdsLayout(const std::string &gdsFile, IndexType cellIdx)
{
    ParserCellGds(_db).parseCellGds(gdsFile, cellIdx);
}

bool IdeaPlaceEx::parseFileBased(int argc, char **argv)
{
    ProgArgs _args = ProgArgsDetails::parseProgArgsCMD(argc, argv);

    // Start message printer timer
    MsgPrinter::startTimer();

    if (!_args.techsimpleFileIsSet())
    {
        ERR("IdeaPlaceEx::%s no techsimple file is given! \n", __FUNCTION__);
        Assert(false);
        return false;
    }
    ParserTechSimple(_db).read(_args.techsimpleFile());
    if (!_args.pinFileIsSet())
    {
        ERR("IdeaPlaceEx::%s no pin file is given! \n", __FUNCTION__);
        Assert(false);
        return false;
    }
    ParserPin(_db).read(_args.pinFile());
    if (!_args.connectionFileIsSet())
    {
        ERR("IdeaPlaceEx::%s no connection is given! \n", __FUNCTION__);
        Assert(false);
        return false;
    }
    ParserConnection(_db).read(_args.connectionFile());
    if (_args.netwgtFileIsSet())
    {
        // If .netwgt file is set, read it
        INF("IdeaPlaceEx::%s Read in the .netwgt ... \n", __FUNCTION__);
        ParserNetwgt(_db).read(_args.netwgtFile());
    }
    else
    {
        // If .netwgt is not set, skip
        INF("IdeaPlaceEx::%s no .netwgt file, skip... \n", __FUNCTION__);
    }

    if (_args.symFileIsSet())
    {
        INF("IdeaPlaceEx:%s Read in .sym ... \n", __FUNCTION__);
        ParserSymFile(_db).read(_args.symFile());
    }

    // Init cells before read in the gds
    if(!_db.initCells())
    {
        ERR("IdeaPlaceEx::%s initializing the cells failed! \n", __FUNCTION__);
        Assert(false);
        return false;
    }
    
    // Parsing the gds files...
    ParserCellGdsDetails::parseAllGdsFiles(_db, _args.gdsFiles());
    return true;
}

bool IdeaPlaceEx::solve()
{
    // Start message printer timer
    MsgPrinter::startTimer();
    // Solve cleaning up tasks for safe...
    for (IndexType cellIdx = 0; cellIdx < _db.numCells(); ++cellIdx)
    {
        _db.cell(cellIdx).calculateCellBBox();
        DBG("cell %d bbox %s \n", cellIdx, _db.cell(cellIdx).cellBBox().toStr().c_str());
    }

    NlpWnconj nlp(_db);
    nlpPtr = &nlp;
    nlp.solve();
#ifdef DEBUG_DRAW
    _db.drawCellBlocks("./debug/after_gr.gds");
#endif //DEBUG_DRAW
    CGLegalizer legalizer(_db);
    legalizer.legalize();
    return true;
}

bool IdeaPlaceEx::outputFileBased(int argc, char **argv)
{
    return true;
}

PROJECT_NAMESPACE_END
