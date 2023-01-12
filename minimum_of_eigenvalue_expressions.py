
# Check whether (sum l^2)^3 / (sum l^3)^2  has minimum when l_1 = positive and all other l's are negative

if __name__ == '__main__':
	
	
	lambda_list =[0,0,0,0]
	for i1 in range(0,10):
		
		for i2 in range(0,100):
			
			for i3 in range(0,100):
				
				for i4 in range(0,100):
				
					lambda_list[0] = (i1+1) * 0.0001
					lambda_list[1] = (i2) * 0.0001
					lambda_list[2] = (i3) * 0.0001
					lambda_list[3] = (i3) * 0.0001
				
				
					lambda_squares =0
					lambda_cubes =0
					for j in range(0, len(lambda_list)):
						lambda_squares = lambda_squares + lambda_list[j]**2
						lambda_cubes = lambda_cubes + lambda_list[j] ** 3
					
					result = lambda_squares**3/ lambda_cubes**2
				
					print(lambda_list, result)
					if result<1 and lambda_list[1]+lambda_list[2]+lambda_list[3]>0.00000000000001:
						input()
					
				