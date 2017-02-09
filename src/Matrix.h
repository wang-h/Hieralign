#ifndef MATRIX_H
#define MATRIX_H

class Matrix
{
	public:
		typedef typename std::vector<double>::reference reference;
		typedef typename std::vector<double>::const_reference const_reference;
		typedef typename std::vector<double>::iterator iterator;
		typedef typename std::vector<double>::const_iterator const_iterator;
		Matrix(): width_(0), height_(0) {};
		
		Matrix(unsigned w, unsigned h) :
			width_(w), height_(h), data_(w*h) {}
		~Matrix(){};
		bool empty() const { return data_.empty(); }
		void resize(unsigned w, unsigned h) {
			data_.resize(w * h);
			width_ = w;
			height_ = h;
		}
		unsigned width() const { return width_; }
		unsigned height() const { return height_; }
		void clear() { data_.clear(); width_=0; height_=0; }
		reference operator()(unsigned i, unsigned j) {
			return data_[offset(i, j)];
		}
		const_reference operator()(unsigned i, unsigned j) const {
			return data_[offset(i, j)];
		}
		iterator begin_col(unsigned j) {
			return data_.begin() + offset(0,j);
		}
		const_iterator begin_col(unsigned j) const {
			return data_.begin() + offset(0,j);
		}
		iterator end_col(unsigned j) {
			return data_.begin() + offset(0,j) + width_;
		}
		const_iterator end_col(unsigned j) const {
			return data_.begin() + offset(0,j) + width_;
		}
		iterator end() { return data_.end(); }
		const_iterator end() const { return data_.end(); }
		void PrintMatrix() {
			for(size_t i = 0; i < width_; i++){ 
				for(size_t j = 0; j < height_; j++){
					printf("%.4f ", data_[offset(i, j)]);
				}
				printf("\n");
			}
		} 
	private:
		inline unsigned offset(unsigned i, unsigned j) const {
			assert(i<width_);
			assert(j<height_);
			return i + j * width_;
		}

		unsigned width_;
		unsigned height_;

		std::vector<double> data_; 
};

#endif // MATRIX_H
