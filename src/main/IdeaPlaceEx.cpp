#include "IdeaPlaceEx.h"
/* Parsers */
#include "parser/ProgArgs.h"
#include "parser/ParserTechSimple.h"
#include "parser/ParserPin.h"
#include "parser/ParserConnection.h"
#include "parser/ParserNetwgt.h"
#include "parser/ParserGds.h"

PROJECT_NAMESPACE_BEGIN

// A global pointer to NLP
// This is for a trick to overly use static members in the NlpWnconj class
NlpWnconj *nlpPtr = nullptr;

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
