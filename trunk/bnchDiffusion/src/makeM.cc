
//build operator matrix
vector<vector<double> > M;

M.resize(n);
for (size_t i=0;i<M.size();++i)
{
  M[i].resize(n,0.0)
}

//functions to fill matrix
for (size_t i=0; i<M.size();++i)
{
  //each i is an equation, operating on phi = 0
}

//computes values along the diagonal
void MD(vector<vector<double> > &M,size_t i)
{
  double val;
  if (i==0)
  {
    val = 0.0; //TODO
  }
  else if (i==M.size()-1)
  {
    val = 0.0; //TODO
  }
  else
  {
    val = 0.0; //TODO
  }
  M[i,i] = val;
}

string whereAmI(const size_t i,const size_t j)
{
  if (i==0)
  {
    if (j==0)
    {
      return string("TopLeft");
    }
    else if (j==J-1)
    {
      return string("BottomLeft");
    }
    else
    {
      return string("MidLeft");
    }
  }
  else if (i==I)
  {
    if (j==0)
    {
      return string("TopRight");
    }
    else if (j==J-1)
    {
      return string("BottomRight");
    }
    else
    {
      return string("MidRight");
    }
  }
  else if (j==0)
  {
    return string("TopMid");
  }
  else if (j==J-1)
  {
    return string("BotMid");
  }
  else
  {
    return string("Interior");
  }
  return string("You should never get this!");
}

//computes values along lower diagonal
void FillM(vector<vector<double> > &M,
             const size_t i, 
             const size_t j, 
             const size_t I, 
             const size_t J,
             const bool offDiag=false)
{
  size_t l = i+j*I;
  double val;
  string loc=whereAmI(i,j);
  switch (loc)
  {
    case string("TopLeft"):
      break;
    case string("TopRight"):
      break;
    case string("BottomLeft"):
      break;
    case string("BottomRight"):
      break;
    case string("MidLeft"):
      break;
    case string("MidRight"):
      break;
    case string("TopMid"):
      break;
    case string("BottomMid"):
      break;
    case string("Interior"):
      break;
    default:
      //should never get here!
  }
}

}
