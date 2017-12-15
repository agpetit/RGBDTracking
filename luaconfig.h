#ifndef LUACONFIG_H
#define LUACONFIG_H

#include <string>

struct lua_State;

class LuaConfig
{
public:
    LuaConfig(const std::string &filename);

    std::string getAsString(const std::string &cmd);
    double getAsNumber(const std::string &cmd);
    std::string returnAsString(const std::string &cmd);
    double returnAsNumber(const std::string &cmd);
private:
    lua_State *L;
};

#endif
