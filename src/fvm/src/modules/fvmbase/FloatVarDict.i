%include <std_map.i>

%include "Field.i"

template<class T>
struct FloatVal
{
private:
  FloatVal();
};

%template(FloatValA) FloatVal<ATYPE_STR>;


template <class T>
class FloatVarDict // : public std::map<string,T>
{
public:

  %extend {
    void setVar(std::string varName, T val)
    {
      FloatVarDict<T>::iterator pos(self->find(varName));
      if (pos != self->end())
        pos->second.constant = val;
      else
        throw CException("uknown var" + varName);
    }

    void __setitem__(std::string varName, T val)
    {
      FloatVarDict<T>::iterator pos(self->find(varName));
      if (pos != self->end())
      {
          pos->second.constant = val;
          pos->second.field = 0;
      }
      else
        throw CException("uknown var" + varName);
    }

    void __setitem__(std::string varName, Field* field)
    {
      FloatVarDict<T>::iterator pos(self->find(varName));
      if (pos != self->end())
        pos->second.field = field;
      else
        throw CException("uknown var" + varName);
    }

    T getVar(std::string varName)
    {
      FloatVarDict<T>::iterator pos(self->find(varName));
      if (pos != self->end())
        return pos->second.constant;
      throw CException("uknown var" + varName);
    }

    FloatVal<T> __getitem__(std::string varName)
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
