#include "lua.hpp"
#include "luaconfig.h"
#include <iostream>
#include <cstdlib>

LuaConfig::LuaConfig(const std::string &filename)
{
    L = lua_open();
    luaL_openlibs(L);
    if (luaL_dofile(L, filename.c_str()))
    {
        std::cerr << "Error loading Lua file '" << filename << "'" << std::endl;
        if (lua_isstring(L, -1))
            std::cerr << lua_tostring(L, -1) << std::endl;
        exit(1);
    }
}

std::string LuaConfig::returnAsString(const std::string &cmd)
{
    if (luaL_dostring(L, cmd.c_str()))
    {
        std::cerr << "Error running Lua code '" << cmd << "'" << std::endl;
        if (lua_isstring(L, -1))
            std::cerr << lua_tostring(L, -1) << std::endl;
        exit(1);
    }
    if (lua_isstring(L, -1))
        return lua_tostring(L, -1);
    return std::string();
}

double LuaConfig::returnAsNumber(const std::string &cmd)
{
    if (luaL_dostring(L, cmd.c_str()))
    {
        std::cerr << "Error running Lua code '" << cmd << "'" << std::endl;
        if (lua_isstring(L, -1))
            std::cerr << lua_tostring(L, -1) << std::endl;
        exit(1);
    }
    if (lua_isnumber(L, -1))
        return lua_tonumber(L, -1);
    return 0.0;
}

std::string LuaConfig::getAsString(const std::string &cmd)
{
    return returnAsString("return " + cmd);
}

double LuaConfig::getAsNumber(const std::string &cmd)
{
    return returnAsNumber("return " + cmd);
}
