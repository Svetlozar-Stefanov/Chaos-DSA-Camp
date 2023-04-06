#include <iostream>
#include <vector>

using namespace std;

float KahanSum(vector<float>& numbers)
{
	float sum = 0;
	float err = 0;
	for (int i = 0; i < numbers.size(); i++)
	{
		float y = numbers[i] - err;
		volatile float t = sum + y;
		volatile float temp = (t - sum);
		err = temp - y;
		sum = t;
	}

	return sum;
}

int main()
{
    vector<float> nums = { 1.5, 2.5, 300.00123, 123434534.4303 };
    std::cout << KahanSum(nums);
}
