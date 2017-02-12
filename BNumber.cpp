#include <vector>
#include <iostream>
#include <string>
#include <cstdio>
#include <cmath>
#include <algorithm>
#include <string>
#include <bits/stdc++.h>
#include <time.h>
using namespace std;
const int base = 1000000000;
class BNumber;
BNumber *mulandsub(BNumber *a,BNumber *b,int q,int j,int n);
vector <int>  add(vector<int>a,vector<int>b)
{
    vector<int>res;
    int carry = 0;
    for (int i = 0; i < max(a.size(),b.size()); ++i) {
        if(a.size() > i)carry += a[i];
        if(b.size() > i)carry += b[i];
        res.push_back(carry%base);
        carry /= base;
    }
    return res;
}
vector <int>  sub(vector<int>a,vector<int>b)
{
    vector<int>res;
    int carry = 0;
    for(int i = 0 ; i <= a.size()-1;i++)
    {
        carry += a[i] - (i < b.size() ? b[i] : 0);
        if (carry < 0) res.push_back(carry+base), carry = -1;
        else res.push_back(carry), carry = 0;
    }
    return res;
}
class BNumber
{
private :
    vector<int>data;
    bool sign;
public:
    /////normalization///////
    void removezeros()
    {
        while (data.size() > 1 && data.back() == 0) data.pop_back();
    }
    //////default constructor////////////
    BNumber()
    {
        sign = false;
    }
    /////////////constructor //////////////
    BNumber(string s)
    {
        if (s.size() == 0)
        {
            data.push_back(0);
        }
        if(s[0]=='-')
        {
            sign = true; // -
            s.erase(0,1);
        }
        else sign = false; // +
        while (s.size()%9 != 0)
            s = '0'+s;
        for (int i=0;i<s.size();i+=9) {
            int v = 0;
            for (int j=i;j<i+9;j++){
                v = v*10+(s[j]-'0');
            }
            data.insert(data.begin(),v);
        }
        this->removezeros();
    }
    ///////////Knuth algorithm ////////////
    //////takes mod parameter , if equal 1 then it is mod else it is div
    BNumber * dimmod(BNumber *a,int mod)
    {
        if(Less(a))
        {
            if(mod)
                return this;
            return new BNumber(0);
        }
        int D = base / (a->data.back() + 1 );
        int n = a->data.size();
        int m = data.size() - n;
        BNumber *d = new BNumber(D);
        BNumber *u = Mul(d);
        if(u->data.size() == m + n)
        {
            u->data.push_back(0);
        }
        BNumber *v = a->Mul(d);
        BNumber *qou = new BNumber();
        qou->data.resize(m+1);
        int indexOfq = m;
        int j = u->data.size()-1;
        for( ;j >= n ; j-- )
        {
            long long q;
            if(u->data[j] == v->data.back())
                q = base-1;
            else
            {
                q = (u->data[j]);
                q *= base;
                q += u->data[j-1];
                q /=v->data.back();
            }
            if(!q)
            {
                indexOfq--;
                continue;
            }
            long long e1 = v->data[n-2];
            e1 *=  q;
            long long e2 = u->data[j];
            e2 *= base;
            e2 += u->data[j-1];
            e2 -=  q * v->data.back();
            e2 *= base;
            e2 += u->data[j-2];
            while( e1 > e2)
            {
                q--;
                long long e3 = v->data[n-2];
                e3 *= q;
                long long e4 = u->data[j];
                e4 *= base;
                e4 += u->data[j-1];
                e4 -=  q * v->data.back();
                e4 *= base;
                e4 += u->data[j-2];
                e1 = e3;
                e2 = e4;
            }
            // * and -
            BNumber *temp = mulandsub(u,v,q,j,n);
            if(temp->getSign())
            {
                q--;
                temp = temp->Add(v);
            }
            int index = 0;
            int t ;
            for( t = j-n ; index  < temp->data.size(); t++)
            {
                u->data[t]=temp->data[index];
                index++;
            }
            for (; index < n+1 ; ++t) {
                u->data[t] = 0;
                index++;
            }
            qou->data[indexOfq] = q;
            indexOfq--;
        }
        if(mod)
        {
            BNumber *r = new BNumber();
            vector<int> o ;
            o.resize(n);
            int index = n-1;
            int count = n;
            for(int q = u->data.size() - m -2 ;count >= 1 ; )
            {
                o[index] = u->data[q];
                index--;
                q--;
                count--;
            }
            r->setData(o);
            r->removezeros();
            r = r->Div(d);
            return r;
        }
        qou->removezeros();
        return qou;
    }
    //////////constructor takes int ////////////
    BNumber(int x) {
        if(x < 0 )
        {
            sign = true;
            x = x * -1;
        }
        if(x == 0)this->data.push_back(0);
        this->setSign(false);
        string s = "";
        while (x > 0) s = char(x%10+'0') + s, x /= 10;
        while (s.size()%9 != 0)
            s = '0'+s;
        for (int i=0;i<s.size();i+=9) {
            int v = 0;
            for (int j=i;j<i+9;j++){
                v = v*10+(s[j]-'0');
            }
            data.insert(data.begin(),v);
        }
        this->removezeros();
    }
    vector<int>getData()
    {
        return data;
    }
    void setData(vector<int> a)
    {
        data = a;
    }

    bool getSign()
    {
        return sign;
    }
    void setSign(bool x)
    {
        sign = x;
    }
    void print()
    {
        if (sign)
        {
            cout<<"-";
        }
        printf("%d", (data.size() == 0) ? 0 : data.back());
        for (int i=data.size()-2;i>=0;i--) printf("%09d", data[i]);
    }
    // get moduls of the BNumber //////////
    BNumber *abs()
    {
        BNumber *abs = new BNumber();
        abs = this;
        abs->setSign(false);
        return abs;
    }
    bool Greater (BNumber *a)
    {
        removezeros();
        a->removezeros();
        if(sign == true && a->getSign() == false)
            return false;
        if(sign == false && a->getSign() == true)
            return true;
        if (a->data.size() > this->data.size())
        {
            return false;
        }else if(a->data.size() < this->data.size())
        {
            return true;
        }else{
            for (int i = a->data.size()-1; i >=0; --i)
            {
                if (a->data[i] > this->data[i])
                {
                    return false;
                }else if(a->data[i] < this->data[i]){
                    return true;
                }
            }
        }
        return false;
    }
    bool Less (BNumber *a)
    {
        removezeros();
        a->removezeros();
        if(sign == true && a->getSign() == false)
            return true;
        if(sign == false && a->getSign() == true)
            return false;
        if (a->data.size() > this->data.size())
        {
            return true;
        }else if(a->data.size() < this->data.size())
        {
            return false;
        }else{
            for (int i = a->data.size()-1; i >=0; --i)
            {
                if (a->data[i] < this->data[i])
                {
                    return false;
                }else if(a->data[i] > this->data[i]){
                    return true;
                }
            }
        }
        return false;
    }
    bool Equal(BNumber *a)
    {
        removezeros();
        a->removezeros();
        if (this->data.size() > a->data.size() || this->data.size() < a->data.size() ||
                sign != a->getSign())
        {
            return false;
        }else{
            for (int i = a->data.size()-1; i >=0; --i)
            {
                if (a->data[i] != this->data[i])
                {
                    return false;
                }
            }
        }
        return true;
    }
     BNumber * Sub (BNumber *a)
    {
        removezeros();
        a->removezeros();
        BNumber *res ;
        a->setSign(!a->getSign());
        res = Add(a);
        return res;

    }
     BNumber *Add (BNumber *a)
    {
        removezeros();
        a->removezeros();
        BNumber *res = new BNumber();
        if(this == 0)return a;
        if(a == 0)return this;
        if(sign == a->getSign())
        {
            res->setSign(sign);
            res->setData(add(a->getData(),data));
        }else{
            if(this->sign == false)
            {
                if(abs()->Greater(a->abs()) || abs()->Equal(a->abs()))
                {
                    res->setSign(false);
                    res->setData(sub(data,a->getData()));
                }else{
                    res->setSign(true);
                    res->setData(sub(a->data,data));
                }
            }else{
                if(abs()->Greater(a->abs()) || abs()->Equal(a->abs()))
                {
                    res->setSign(true);
                    res->setData(sub(data,a->getData()));
                }else{
                    res->setSign(false);
                    res->setData(sub(a->data,data));
                }
            }
        }
        res->removezeros();
        return res;
    }
    BNumber *Mul(BNumber *a)
    {
        removezeros();
        a->removezeros();
        BNumber *ans = new BNumber();
        if(this->getSign() != a->getSign())
        {
            ans->setSign(true);
        }else{
            ans->setSign(false);
        }

        unsigned long long  carry = 0;
        ans->data.assign(a->data.size()+data.size(), 0);
        for(int i=0;i<=data.size()-1;i++){
            carry = 0ll;
            for (int j=0;j<a->data.size() || carry > 0;j++) {
                unsigned long long s = ans->data[i+j] + carry + (unsigned long long)data[i]*(j < a->data.size() ? (unsigned long long)a->data[j] : 0ll);
                ans->data[i+j] = s % base;
                carry = s / base;
            }
        }
        ans->removezeros();
        return ans;
    }
    BNumber * Div(BNumber *a) {
        if(a->data.size() > 1)
        {
            return dimmod(a,0);
        }
        long long b = a->data[0];
        BNumber *v = new BNumber();
        v->data.resize(data.size());
        long long cur = 0;
        for(int i = this->data.size()-1;i >= 0 ; i--)
        {
            cur = (cur % b) * base + data[i];
            v->data[i] = cur / b;
        }
        v->removezeros();
        return v;
    }
     BNumber  *Mod (BNumber *a) {
        if(a->data.size() > 1)
        {
            return this->dimmod(a,1);
        }
        long long b = a->data[0];
        BNumber *v = new BNumber ();
        long long cur = 0;
        for (int i = data.size()-1; i >= 0 ; --i) {
            cur *= base;
            cur += data[i];
            cur %= b;
        }
        v->data.push_back(cur);
        return v;

    }
    string isprime(int iteration)
    {
        int i;
        BNumber *ONE = new BNumber(1);
        BNumber *TWO = new BNumber(2);
        BNumber *ZERO = new BNumber(0);
        BNumber *k = Sub(ONE);
        BNumber *s = k;
        if (Less(TWO))
        {
            return "No";
        }
        if (!Equal(TWO) && Mod(TWO)->Equal(ZERO))
        {
            return "No";
        }
        while (s->Mod(TWO)->Equal(ZERO))
        {
            s = s->Div(TWO);
        }
        for (i = 0; i < iteration; i++)
        {
            BNumber *a= new BNumber(rand());
            a = a->Mod(k);
            a = a->Add(ONE);
            BNumber *temp = s;
            BNumber * mod = a->power(temp,this);
            while (!temp->Equal(k) && !mod->Equal(ONE) && !mod->Equal(k))
            {
                mod = mod->mulmod(mod,this);

                temp = temp->Mul(TWO);
            }
            if (!mod->Equal(k) && temp->Div(TWO)->Equal(ZERO))
            {
                return "No";
            }
        }
        return "Yes";
    }
    BNumber *mulmod(BNumber *a, BNumber *mod)
    {
        a = a->Mod(mod);
        BNumber * out = Mod(mod);
        out = out->Mul(a);
        return out->Mod(mod);
    }
    BNumber * extendedEucluid(BNumber *m)
    {
        BNumber *ONE = new BNumber(1);
        BNumber *TWO = new BNumber(2);
        BNumber *ZERO = new BNumber(0);
        BNumber *A1 = new BNumber(1);
        BNumber *A2 = new BNumber(0);
        BNumber *A3 = m;
        BNumber *B1 = new BNumber(0);
        BNumber *B2 = new BNumber(1);
        BNumber *B3 = this;
L: if(B3->Equal(ZERO))return this->gcd(m);
        if(B3->Equal(ONE))
        {
            if(B2->getSign() == true)
            {
                B2 = B2->Add(m);
            }
            return B2;
        }
        BNumber *Q = A3->Div(B3);
        BNumber *T1 = A1->Sub(Q->Mul(B1));
        BNumber *T2 = A2->Sub(Q->Mul(B2));
        BNumber  *T3= A3->Sub(Q->Mul(B3));
        A1 =  B1;
        A2 = B2 ;
        A3 = B3;
        B1 = T1;
        B2 = T2 ;
        B3 = T3 ;
        goto L;
    }
    BNumber *gcd(BNumber *b) {
        BNumber *ZERO = new BNumber(0);
        BNumber *a = this;
        while (b->Greater(ZERO)) {
            BNumber *r = a->Mod(b);
            a = b;
            b = r;
        }
        return a;
    }
    BNumber *power(BNumber *a ,BNumber *mod)
    {
        BNumber *ZERO = new BNumber(0);
        BNumber *ONE = new BNumber(1);
        BNumber *TWO = new BNumber(2);
        BNumber *x = new BNumber(1);
        BNumber *y = this;
        while (a->Greater(ZERO))
        {
            if (a->Mod(TWO)->Equal(ONE))
            {
                x = x->Mul(y)->Mod(mod);
            }
            y =y->Mul(y)->Mod(mod);
            a = a->Div(TWO);
        }
        return x->Mod(mod);
    }
};
BNumber *phiofn(BNumber *p,BNumber *q)
{
    BNumber *ONE = new BNumber(1);
    return q->Sub(ONE)->Mul(p->Sub(ONE));
}
void program()
{
    string input1;
    cin >> input1;
    string input2;
    cin >> input2;
    string input3;
    cin >> input3;
    string p1 = input1.erase(0,2);
    string q1 = input2.erase(0,2);
    string e1 = input3.erase(0,2);
    BNumber *p = new BNumber(p1);
    BNumber *q = new BNumber(q1);
    BNumber *e = new BNumber(e1);
    BNumber *pn = phiofn(p,q);
    BNumber *n = p->Mul(q);
    string operation;
    cin >> operation;
    while (operation != "Quit") {
        if(operation == "IsPPrime")
        {
            cout<< p->isprime(1) <<endl;
        }else if(operation =="IsQPrime")
        {
            cout << q->isprime(1) << endl;
        }else if(operation == "PrintN")
        {
            n->print();cout<<endl;
        }else if(operation =="PrintD")
        {
            BNumber *d = e->extendedEucluid(pn);
            d->print();cout<<endl;
        }else if(operation == "PrintPhi")
        {
            pn->print();cout<<endl;
        }else if(operation.substr(0,13) == "EncryptPublic")
        {
            string m = operation.substr(15,-1);
            m.pop_back();
            BNumber *message = new BNumber(m);
            BNumber *out = message->power(e,n);
            out->print();cout<<endl;
        }else if(operation.substr(0,14) == "EncryptPrivate")
        {
            BNumber *d = e->extendedEucluid(pn);
            string m = operation.substr(16,-1);
            m.pop_back();
            BNumber *message = new BNumber(m);
            BNumber *out =  message->power(d,n);
            out->print();cout<<endl;
        }else if(operation == "Quit")
        {
            break;
        }
        cin >> operation;
    }

}
BNumber * mulandsub(BNumber *a,BNumber *b,int z,int j,int n)
{
    BNumber *q = new BNumber(z);
    BNumber * temp = b->Mul(q);
    BNumber *temp1=new BNumber();
    vector<int>subu;
    subu.resize(n+1);
    int iu = n;
    for(int o = j ; o >= j - n ; o--)
    {
        subu[iu] = a->getData()[o];
        iu--;
    }
    temp1->setData(subu);
    BNumber * out = new BNumber();
    out = temp1->Sub(temp);
    return out;
}
int main()
{
    freopen("i","r",stdin);

    program();
    program();
    return 0;
}
