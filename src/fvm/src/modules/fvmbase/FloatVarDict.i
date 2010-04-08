%include <std_map.i>


template<class T>
struct FloatVal
{
private:
  FloatVal();
};

%template(FloatValA) FloatVal< ATYPE_STR >;

%include "atype.i"

template <class T>
class FloatVarDict // : public std::map<string,T>
{
public:

  %extend {
    void setVar(std::string varName, T val)
    {
      FloatVarDict< T >::iterator pos(self->find(varName));
      if (pos != self->end())
        pos->second.constant = val;
      else
        throw CException("uknown var" + varName);
    }

    std::vector<std::string>
      getKeys()
    {
        std::vector<string> keys;
        for(FloatVarDict< T >::iterator pos(self->begin());
            pos != self->end();
            ++pos)
        {
            keys.push_back(pos->first);
        }
        return keys;
    }
    
    void __setitem__(std::string varName, T val)
    {
      FloatVarDict< T >::iterator pos(self->find(varName));
      if (pos != self->end())
      {
          pos->second.constant = val;
          pos->second.field = 0;
      }
      else
        throw CException("uknown var" + varName);
    }

#ifdef USING_ATYPE_TANGENT
    void __setitem__(std::string varName, double val)
    {
      FloatVarDict< T >::iterator pos(self->find(varName));
      if (pos != self->end())
      {
          pos->second.constant = T(val);
          pos->second.field = 0;
      }
      else
        throw CException("uknown var" + varName);
    }
#endif
    
#ifdef USING_ATYPE_PC
    void __setitem__(std::string varName, double val)
    {
      FloatVarDict< T >::iterator pos(self->find(varName));
      if (pos != self->end())
      {
          pos->second.constant = T(val);
          pos->second.field = 0;
      }
      else
        throw CException("uknown var" + varName);
    }
#endif
    
    void __setitem__(std::string varName, Field* field)
    {
      FloatVarDict< T >::iterator pos(self->find(varName));
      if (pos != self->end())
        pos->second.field = field;
      else
        throw CException("uknown var" + varName);
    }

    T getVar(std::string varName)
    {
      FloatVarDict< T >::iterator pos(self->find(varName));
      if (pos != self->end())
        return pos->second.constant;
      throw CException("uknown var" + varName);
    }

    FloatVal< T > __getitem__(std::string varName)
    {
      FloatVarDict< T >::iterator pos(self->find(varName));
      if (pos != self->end())
        return pos->second;
      throw CException("uknown var" + varName);
    }
  }

};

//%template(mapA) std::map<string,ATYPE_STR>;
%template(FloatVarDictA) FloatVarDict< ATYPE_STR >;
