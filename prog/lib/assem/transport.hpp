#ifndef BF_TRANSPORT_HPP
#define BF_TRANSPORT_HPP

#include <vector>
#include "grid.hpp"

namespace bf{

class ITransport{
public:
	virtual void set_velocity(const std::vector<double>& vel) = 0;
	virtual std::vector<double> compute_fluxes(const std::vector<double>& u) const = 0;
};

class UpwindTransport: public ITransport{
public:
	void set_velocity(const std::vector<double>& vel) override;
	std::vector<double> compute_fluxes(const std::vector<double>& u) const override;
private:
	std::vector<double> _velocity;
};

class TvdSuperbeeTransport: public ITransport{
public:
	void set_velocity(const std::vector<double>& vel) override;
	std::vector<double> compute_fluxes(const std::vector<double>& u) const override;
private:
	std::vector<double> _velocity;
};

class TvdMcTransport: public ITransport{
public:
	void set_velocity(const std::vector<double>& vel) override;
	std::vector<double> compute_fluxes(const std::vector<double>& u) const override;
private:
	std::vector<double> _velocity;
};

class TvdVanLeerTransport: public ITransport{
public:
	void set_velocity(const std::vector<double>& vel) override;
	std::vector<double> compute_fluxes(const std::vector<double>& u) const override;
private:
	std::vector<double> _velocity;
};

class TvdMinModTransport: public ITransport{
public:
	void set_velocity(const std::vector<double>& vel) override;
	std::vector<double> compute_fluxes(const std::vector<double>& u) const override;
private:
	std::vector<double> _velocity;
};

}
#endif
