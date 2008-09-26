%include <std_map.i>

template <class T>
class FloatVarDict // : public std::map<string,T>
{
public:

  %extend {
    void setVar(std::string varName, T val)
    {
      FloatVarDict<T>::iterator pos(self->find(varName));
      if (pos != self->end())
        pos->second = val;
      else
        throw CException("uknown var" + varName);
    }

    T getVar(std::string varName)
    {
      FloatVarDict<T>::iterator pos(self->find(varName));
      if (pos != self->end())
        return pos->second;
      throw CException("uknown var" + varName);
    }
  }

};

//%template(mapA) std::map<string,ATYPE_STR>;
%template(FloatVarDictA) FloatVarDict<ATYPE_STR>;
