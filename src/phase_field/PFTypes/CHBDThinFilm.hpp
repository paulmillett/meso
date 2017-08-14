
# ifndef CHBDTHINFILM_H
# define CHBDTHINFILM_H

# include "CHBD.hpp"

class CHBDThinFilm: public CHBD {

    public:

        CHBDThinFilm(const CommonParams&, const GetPot&);
        ~CHBDThinFilm();
        virtual void initPhaseField(); // override
        virtual void updatePhaseField(); // override

    protected:
        void makeWalls();
        int thickness;

};

# endif  // CHBDTHINFILM_H
