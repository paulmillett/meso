
# ifndef CHAPP_H
# define CHAPP_H

# include <vector>
# include "../../base/MesoBase.hpp"
# include "Params.hpp"
# include "CHSystem.hpp"

class CHApp : public MesoBase {

private:

	int current_step;
	Params p;
	CHSystem* ch_system;
	
public:

	CHApp(const GetPot&);
	~CHApp();
	void initSystem();
	void stepForward(int);
	void writeOutput(int);

};

# endif  // CHAPP_H
