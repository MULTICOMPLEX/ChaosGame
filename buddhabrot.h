
class Buddhabrot 
{
	public:
	Buddhabrot(int height, int width) : m_height(height), m_width(width) {
		image.resize((m_width + 1) * (m_height + 1) * 4);
	}
	void gen_fractal();
	std::vector<unsigned char> image; 

	unsigned int get_width() const;
	unsigned int get_height() const;
	
protected:
	unsigned int m_height;
	unsigned int m_width;
};

