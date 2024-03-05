def out_of_bounds(x,y, angle_in):
	x_min = -10
	x_max = 10
	y_min = -10
	y_max = 10

	if x < x_min:
		y -= (x - x_min)*math.tan(angle_in)
		print(x_min, y, "x_min")
		return True, x_min, y, 1
	elif x > x_max:
		y -= (x - x_max)*math.tan(angle_in)
		print(x_max, y, "x_max")
		return True, x_max, y, 1
	elif y < y_min:
		x -= (y - y_min)/math.tan(angle_in)
		print(x, y_min, "y_min")
		return True, x, y_min, 0
	elif y > y_max:
		x -= (y - y_max)/math.tan(angle_in)
		print(x, y_max, "y_max")
		return True, x, y_max, 0
	else:
		return False, x, y,0

def test():
	x2,y2 = 9,11
	x1,y1 = 9,9
	print("Out of bound", x2, y2)
	dx = x2-x1
	dy = y2-y1
	angle = math.atan2(dy,dx)
	dr = math.sqrt(dx**2+dy**2)
	OoB, x_wall, y_wall, reflect_corr = out_of_bounds(x2,y2,angle)
	dr_wall = math.sqrt((x2-x_wall)**2 + (y2-y_wall)**2)
	new_angle = math.pi*reflect_corr - angle
	x = x_wall + (dr-dr_wall)*math.cos(new_angle)
	y = y_wall + (dr-dr_wall)*math.sin(new_angle)
	print("New", x, y)

# test()