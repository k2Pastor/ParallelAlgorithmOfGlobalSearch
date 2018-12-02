#define _USE_MATH_DEFINES
#include <iostream>
#include <mutex>
#include <thread>
#include <queue>
#include <list>
#include <iterator>
#include <cmath>
#include <ctime>

using namespace std;

mutex m_lock;
mutex tp_lock;
mutex queue_lock;

struct TestPoint
{
	double x;
	double z;
	TestPoint(double _x = 0, double _z = 0) : x(_x), z(_z) { }
};

struct Interval
{
	double R1;
	TestPoint* pLeft;
	TestPoint* pRight;
	Interval(double _R = 0, TestPoint *tpLeft = NULL, TestPoint *tpRight = NULL)
		: R1(_R), pLeft(tpLeft), pRight(tpRight) { }
};

bool operator<(const Interval& i1, const Interval& i2) { return (i1.R1 < i2.R1) ? true : false; }

double f(double x)
{
return 2 * (x - 3.0) * (x - 3.0) + exp(x * x / 2.0); // совпадение;
//return sin(x) + sin((10 * x) / 3); // совпадение;
//return x * x;
//return ((3.0 * x - 1.4) * sin(18.0 * x)); // совпадение;
//return -1.0 * (x + sin(x)) * exp(-1.0 * x * x); // совпадение;
//return sin(x) + sin((10 * x) / 3) + log(x) - 0.84 * x + 3; // совпадение;
//return -1.0 * sin(2 * M_PI * x) * exp(-1.0 * x); // почти совпадение;
//return (x * x - 5 * x + 6) / (x * x + 1); // совпадение;
//return -x + sin(3 * x) - 1; // совпадение;
} 

/*double f(double x)
{
	double sum = 0.0;
	for (int k = 1; k <= 5; k++)
	{
		sum += k * sin((k + 1) * x + k);
	}

	return -1.0 * sum;
} */

double ComputeR(const TestPoint& tpLeft, const TestPoint& tpRight, double _m)
{
	double diffX = tpRight.x - tpLeft.x;
	double diffZ = tpRight.z - tpLeft.z;
	return _m * diffX + diffZ * diffZ / (_m * diffX) - 2 * (tpRight.z + tpLeft.z);
}

TestPoint* InsertUp_List(list<TestPoint> &ltp, TestPoint &tpk)
{
	list<TestPoint>::iterator itLeft, itRight;
	itLeft = itRight = ltp.begin();

	while ((itRight != ltp.end()) && (itRight->x < tpk.x))
	{
		itLeft = itRight;
		itRight++;
	}

	ltp.insert(itRight, tpk);
	itLeft++;

	return &(*itLeft);
}

void thread_tpk(Interval &i, TestPoint &tp, double m)
{
	double tpk = 0.5 * (i.pRight->x + i.pLeft->x) - ((i.pRight->z - i.pLeft->z) / (2.0 * m));
	tp.x = tpk;
	tp.z = f(tpk);
}

int main()
{
	list<TestPoint> testPoints; // точки испытаний;
	list<TestPoint>::iterator itLeft, itRight;
	priority_queue<Interval> Queue;

	TestPoint tp1, tp2, tp3, tp4, tp5, tp6, tpk, *tpt, DotOfGM;
	Interval CharacteristicInterval1, CharacteristicInterval2, CharacteristicInterval3, CharacteristicInterval4;
	double accuracy; // точность алгоритма;
	int iterations; // количество итераций;
	int k = 0; // количество испытаний;
	double M;
	double m = -1.0; // константа Липшица;
	double r = 2.0; // заданный параметр метода;
	double GlobalMin, DotOfGlobalMin;
	double timeStart, timeEnd;


	cout << "Enter the ends of the segment [a;b] : " << endl;
	cin >> tp1.x >> tp2.x;
	cout << "Enter number of the iterations: " << endl;
	cin >> iterations;
	cout << "Enter value of the accuracy: " << endl;
	cin >> accuracy;

	timeStart = clock();

	tp1.z = f(tp1.x);
	tp2.z = f(tp2.x);

	testPoints.push_back(tp1);
	testPoints.push_back(tp2);

	k = 1;

	do
	{
		M = 0.0;
		itLeft = itRight = testPoints.begin();
		++itRight;

		double Old_m = m;

		while (itRight != testPoints.end())
		{
			double max_M = fabs((itRight->z - itLeft->z) / (itRight->x - itLeft->x));
			if (M < max_M)
			{
				M = max_M;
			}

			itRight++;
			itLeft++;
		}

		if (M > 0.0)
		{
			m = r * M;
		}
		else {
			m = 1.0;
		}




		if (Old_m != m)
		{
			Queue = priority_queue<Interval>();

			itLeft = itRight = testPoints.begin();
			itRight++;

			while (itRight != testPoints.end())
			{
				Queue.push(Interval(ComputeR(*itLeft, *itRight, m), &(*itLeft), &(*itRight)));
				++itLeft;
				++itRight;
			}

		}

		while (itRight != testPoints.end())
		{
			Queue.push(Interval(ComputeR(*itLeft, *itRight, m), &(*itLeft), &(*itRight)));
			++itLeft;
			++itRight;
		}

		if (itRight == testPoints.end())
		{
			(itRight--);
		}

	switch (k)
	{

		case 1: 
		{
			CharacteristicInterval1 = Queue.top();
			Queue.pop();

			thread thr1(thread_tpk, ref(CharacteristicInterval1), ref(tp3), m);

			thr1.join();

			tpt = InsertUp_List(testPoints, tp3);
			Queue.push(Interval(ComputeR(*CharacteristicInterval1.pLeft, *tpt, m), CharacteristicInterval1.pLeft, tpt));
			Queue.push(Interval(ComputeR(*tpt, *CharacteristicInterval1.pRight, m), tpt, CharacteristicInterval1.pRight));

			k++;

			break;
		}

		case 2:
		{
			CharacteristicInterval1 = Queue.top();
			Queue.pop();
			CharacteristicInterval2 = Queue.top();
			Queue.pop();

			thread thr1(thread_tpk, ref(CharacteristicInterval1), ref(tp3), m);
			thread thr2(thread_tpk, ref(CharacteristicInterval2), ref(tp4), m);

			thr1.join();
			thr2.join();

			tpt = InsertUp_List(testPoints, tp3);
			Queue.push(Interval(ComputeR(*CharacteristicInterval1.pLeft, *tpt, m), CharacteristicInterval1.pLeft, tpt));
			Queue.push(Interval(ComputeR(*tpt, *CharacteristicInterval1.pRight, m), tpt, CharacteristicInterval1.pRight));

			tpt = InsertUp_List(testPoints, tp4);
			Queue.push(Interval(ComputeR(*CharacteristicInterval2.pLeft, *tpt, m), CharacteristicInterval2.pLeft, tpt));
			Queue.push(Interval(ComputeR(*tpt, *CharacteristicInterval2.pRight, m), tpt, CharacteristicInterval2.pRight));

			k += 2;

			break;
		}

		case 3:
		{
			CharacteristicInterval1 = Queue.top();
			Queue.pop();
			CharacteristicInterval2 = Queue.top();
			Queue.pop();
			CharacteristicInterval3 = Queue.top();
			Queue.pop();

			thread thr1(thread_tpk, ref(CharacteristicInterval1), ref(tp3), m);
			thread thr2(thread_tpk, ref(CharacteristicInterval2), ref(tp4), m);
			thread thr3(thread_tpk, ref(CharacteristicInterval3), ref(tp5), m);

			thr1.join();
			thr2.join();
			thr3.join();

			tpt = InsertUp_List(testPoints, tp3);
			Queue.push(Interval(ComputeR(*CharacteristicInterval1.pLeft, *tpt, m), CharacteristicInterval1.pLeft, tpt));
			Queue.push(Interval(ComputeR(*tpt, *CharacteristicInterval1.pRight, m), tpt, CharacteristicInterval1.pRight));

			tpt = InsertUp_List(testPoints, tp4);
			Queue.push(Interval(ComputeR(*CharacteristicInterval2.pLeft, *tpt, m), CharacteristicInterval2.pLeft, tpt));
			Queue.push(Interval(ComputeR(*tpt, *CharacteristicInterval2.pRight, m), tpt, CharacteristicInterval2.pRight));

			tpt = InsertUp_List(testPoints, tp5);
			Queue.push(Interval(ComputeR(*CharacteristicInterval3.pLeft, *tpt, m), CharacteristicInterval3.pLeft, tpt));
			Queue.push(Interval(ComputeR(*tpt, *CharacteristicInterval3.pRight, m), tpt, CharacteristicInterval3.pRight));

			k += 3;

			break;
		}

		default:
		{
			CharacteristicInterval1 = Queue.top();
			Queue.pop();
			CharacteristicInterval2 = Queue.top();
			Queue.pop();
			CharacteristicInterval3 = Queue.top();
			Queue.pop();
			CharacteristicInterval4 = Queue.top();
			Queue.pop();

			thread thr1(thread_tpk, ref(CharacteristicInterval1), ref(tp3), m);
			thread thr2(thread_tpk, ref(CharacteristicInterval2), ref(tp4), m);
			thread thr3(thread_tpk, ref(CharacteristicInterval3), ref(tp5), m);
			thread thr4(thread_tpk, ref(CharacteristicInterval4), ref(tp6), m);

			thr1.join();
			thr2.join();
			thr3.join();
			thr4.join();

			tpt = InsertUp_List(testPoints, tp3);
			Queue.push(Interval(ComputeR(*CharacteristicInterval1.pLeft, *tpt, m), CharacteristicInterval1.pLeft, tpt));
			Queue.push(Interval(ComputeR(*tpt, *CharacteristicInterval1.pRight, m), tpt, CharacteristicInterval1.pRight));

			tpt = InsertUp_List(testPoints, tp4);
			Queue.push(Interval(ComputeR(*CharacteristicInterval2.pLeft, *tpt, m), CharacteristicInterval2.pLeft, tpt));
			Queue.push(Interval(ComputeR(*tpt, *CharacteristicInterval2.pRight, m), tpt, CharacteristicInterval2.pRight));

			tpt = InsertUp_List(testPoints, tp5);
			Queue.push(Interval(ComputeR(*CharacteristicInterval3.pLeft, *tpt, m), CharacteristicInterval3.pLeft, tpt));
			Queue.push(Interval(ComputeR(*tpt, *CharacteristicInterval3.pRight, m), tpt, CharacteristicInterval3.pRight));

			tpt = InsertUp_List(testPoints, tp6);
			Queue.push(Interval(ComputeR(*CharacteristicInterval4.pLeft, *tpt, m), CharacteristicInterval4.pLeft, tpt));
			Queue.push(Interval(ComputeR(*tpt, *CharacteristicInterval4.pRight, m), tpt, CharacteristicInterval4.pRight));

			k += 4;

			break;

		}
	}

	} while ((fabs(CharacteristicInterval1.pRight->x - CharacteristicInterval1.pLeft->x) > accuracy) || (fabs(CharacteristicInterval2.pRight->x - CharacteristicInterval2.pLeft->x) > accuracy) || (fabs(CharacteristicInterval3.pRight->x - CharacteristicInterval3.pLeft->x) > accuracy) || (fabs(CharacteristicInterval4.pRight->x - CharacteristicInterval4.pLeft->x) > accuracy) && ((k - 2) < iterations));

	itLeft = testPoints.begin();
	GlobalMin = itLeft->z;
	DotOfGlobalMin = itLeft->x;
	itLeft++;

	while (itLeft != testPoints.end())
	{
		if (GlobalMin > itLeft->z)
		{
			GlobalMin = itLeft->z;
			DotOfGlobalMin = itLeft->x;
		}

		itLeft++;
	}

	timeEnd = clock();

	cout.precision(9);
	cout << " Global minimum is: " << std::fixed << GlobalMin << endl;
	cout << " Dot of global minimum is: " << DotOfGlobalMin << endl;
	cout << " Number of experiments is: " << k << endl;
	cout << "Time of algorithm is: " << (timeEnd - timeStart) / CLOCKS_PER_SEC << endl;

	return 0;
}