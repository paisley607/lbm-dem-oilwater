/* =====================================================================================
   parameterConfig.h —— C++ 侧的参数文件读取器（与 config/load_parameters.sh 规则一致）

   【它解决什么问题】
     run_production.sh 已经用 bash 读过一遍 .env 参数，但主程序 oilwaterDrive 还需要
     两相、润湿、入口控制、驱动协议等几十个不走命令行的参数。为了不把它们一个个塞进
     命令行，主程序自己按同一套规则再读一遍 .env 文件。

   【读取规则】（和 bash 版完全一样，三条）
     1. PARAMETER_DEFAULTS_FILE（默认 config/defaults.env）是合法参数名清单 + 默认值；
     2. CONFIG_FILE 指定的实验文件可以用 BASE_CONFIG 继承别的文件，先父后子，后者覆盖前者；
        出现 defaults 里没有的参数名、或循环继承，直接抛异常中止；
     3. 同名环境变量优先级最高（text() 每次先查 getenv）。

   【怎么用】
     water_fangsha::ParameterConfig cfg;          // 构造时读完全部文件
     double tau = cfg.real("TAU_W_TARGET");       // 实数
     long   n   = cfg.integer("RESTART_STAGES");  // 整数
     std::string mode = cfg.text("INLET_CONTROL");// 原始文本
     std::vector<double> v = cfg.realList("X");   // 逗号分隔的实数列表
     任何格式错误都抛 std::runtime_error，主程序在 main 里捕获并打印后退出。
   ===================================================================================== */
#ifndef WATER_FANGSHA_PARAMETER_CONFIG_H
#define WATER_FANGSHA_PARAMETER_CONFIG_H

#include <cerrno>     // errno / ERANGE：strtod、strtol 溢出检测
#include <cmath>      // std::isfinite
#include <cstdlib>    // std::getenv、std::strtod、std::strtol、realpath
#include <fstream>    // std::ifstream 逐行读文件
#include <limits.h>   // PATH_MAX
#include <map>
#include <sstream>    // 解析逗号分隔列表
#include <stdexcept>  // std::runtime_error
#include <string>
#include <vector>

namespace water_fangsha {

/// 参数表：构造时把 defaults + 实验文件（含 BASE_CONFIG 链）合并进 values，
/// 之后按名字取值，取值时环境变量优先。
class ParameterConfig {
public:
    /// 构造：读 PARAMETER_DEFAULTS_FILE（登记合法键），再读 CONFIG_FILE（覆盖）。
    /// CONFIG_FILE 未设置、为空或就是 defaults 本身时，只读 defaults。
    ParameterConfig()
    {
        char const* defaults = std::getenv("PARAMETER_DEFAULTS_FILE");
        char const* profile = std::getenv("CONFIG_FILE");
        std::string defaultsPath = defaults && *defaults ? defaults : "config/defaults.env";
        std::map<std::string,bool> loading;   // 正在递归加载中的文件集合，用来检测 BASE_CONFIG 循环
        loadFile(defaultsPath, loading, true);
        if (profile && *profile && canonicalPath(profile) != canonicalPath(defaultsPath)) {
            loadFile(profile, loading, false);
        }
    }

    /// 取原始文本。优先级：非空环境变量 > 文件值；都没有则抛 "missing parameter"。
    std::string text(std::string const& name) const
    {
        char const* environment = std::getenv(name.c_str());
        if (environment && *environment) return std::string(environment);
        std::map<std::string,std::string>::const_iterator found = values.find(name);
        if (found == values.end()) {
            throw std::runtime_error("missing parameter " + name);
        }
        return found->second;
    }

    /// 取有限实数。拒绝：空串、尾随非数字字符（如 "0.8abc"）、溢出、inf/nan。
    double real(std::string const& name) const
    {
        std::string raw = text(name);
        char* end = 0;
        errno = 0;
        double value = std::strtod(raw.c_str(), &end);
        if (end == raw.c_str() || *end != '\0' || errno == ERANGE || !std::isfinite(value)) {
            throw std::runtime_error("invalid real parameter " + name + "=" + raw);
        }
        return value;
    }

    /// 取十进制整数。拒绝：空串、尾随字符（含 "5.0"）、溢出。
    long integer(std::string const& name) const
    {
        std::string raw = text(name);
        char* end = 0;
        errno = 0;
        long value = std::strtol(raw.c_str(), &end, 10);
        if (end == raw.c_str() || *end != '\0' || errno == ERANGE) {
            throw std::runtime_error("invalid integer parameter " + name + "=" + raw);
        }
        return value;
    }

    /// 取逗号分隔的实数列表（如 "1.0,1.5,2.0"）。每一项按 real() 的规则检查；空列表报错。
    std::vector<double> realList(std::string const& name) const
    {
        std::string raw = text(name);
        std::stringstream stream(raw);
        std::vector<double> valuesList;
        std::string item;
        while (std::getline(stream, item, ',')) {
            char* end = 0;
            errno = 0;
            double value = std::strtod(item.c_str(), &end);
            if (end == item.c_str() || *end != '\0' || errno == ERANGE || !std::isfinite(value)) {
                throw std::runtime_error("invalid list parameter " + name + "=" + raw);
            }
            valuesList.push_back(value);
        }
        if (valuesList.empty()) throw std::runtime_error("empty parameter " + name);
        return valuesList;
    }

    /// 返回"最终生效值"的完整快照（文件值 + 环境变量覆盖），主程序用它把参数写进输出目录存档。
    std::map<std::string,std::string> snapshot() const
    {
        std::map<std::string,std::string> result(values);
        for (std::map<std::string,std::string>::const_iterator it=values.begin(); it!=values.end(); ++it) {
            char const* environment=std::getenv(it->first.c_str());
            if (environment && *environment) result[it->first]=environment;
        }
        return result;
    }

private:
    /// 去掉两端空白和成对引号（"abc" / 'abc' → abc），与 bash 版 trim_parameter_value 一致。
    static std::string trim(std::string value)
    {
        std::string::size_type first = value.find_first_not_of(" \t\r\n");
        if (first == std::string::npos) return std::string();
        std::string::size_type last = value.find_last_not_of(" \t\r\n");
        value = value.substr(first, last-first+1);
        if (value.size() >= 2 && ((value[0] == '\'' && value[value.size()-1] == '\'')
                                  || (value[0] == '"' && value[value.size()-1] == '"'))) {
            value = value.substr(1, value.size()-2);
        }
        return value;
    }

    /// 把路径规范化（解析符号链接、../ 等），用来判断"两个路径是否同一个文件"。
    /// 路径不存在时原样返回，留给后面的打开失败给出可读的错误信息。
    static std::string canonicalPath(std::string const& path)
    {
        char resolved[PATH_MAX];
        char const* result=realpath(path.c_str(), resolved);
        return result ? std::string(result) : path;
    }

    /// 取路径的目录部分（没有分隔符时返回 "."）。
    static std::string directoryName(std::string const& path)
    {
        std::string::size_type slash=path.find_last_of("/\\");
        return slash == std::string::npos ? std::string(".") : path.substr(0,slash);
    }

    /// 解析 BASE_CONFIG：绝对路径原样用；相对路径相对于"声明它的那个文件"所在目录。
    static std::string resolveBase(std::string const& file, std::string base)
    {
        if (base.empty() || base[0]=='/') return base;
        return directoryName(file) + "/" + base;
    }

    /// 读一个参数文件并递归处理它的 BASE_CONFIG。
    ///   registerKeys=true  ：这是 defaults 文件，遇到的每个 KEY 都登记为合法键；
    ///   registerKeys=false ：这是实验/基础文件，KEY 必须已登记，否则视为拼写错误。
    /// 覆盖顺序：先递归加载父文件写入 values，再把本文件的值写入（本文件覆盖父文件）。
    void loadFile(std::string const& path, std::map<std::string,bool>& loading, bool registerKeys)
    {
        std::string canonical=canonicalPath(path);
        if (loading.find(canonical)!=loading.end())
            throw std::runtime_error("cyclic BASE_CONFIG inheritance at " + canonical);
        std::ifstream input(canonical.c_str());
        if (!input) throw std::runtime_error("cannot open parameter file " + path);
        loading[canonical]=true;
        std::map<std::string,std::string> local;   // 本文件自己的 KEY=VALUE，先暂存，等父文件加载完再合并
        std::string line;
        while (std::getline(input, line)) {
            line = trim(line);
            if (line.empty() || line[0] == '#') continue;
            std::string::size_type equals = line.find('=');
            if (equals == std::string::npos) throw std::runtime_error("invalid parameter line: " + line);
            std::string name = trim(line.substr(0, equals));
            std::string value = trim(line.substr(equals+1));
            if (name.empty()) throw std::runtime_error("empty parameter name in " + path);
            if (!registerKeys && allowedKeys.find(name)==allowedKeys.end())
                throw std::runtime_error("unknown parameter " + name + " in " + path);
            if (registerKeys) allowedKeys[name]=true;
            local[name] = value;                   // 同一文件里重复出现时，最后一次赋值生效
        }
        std::map<std::string,std::string>::const_iterator base=local.find("BASE_CONFIG");
        if (base!=local.end() && !base->second.empty())
            loadFile(resolveBase(canonical,base->second),loading,false);
        for (std::map<std::string,std::string>::const_iterator it=local.begin(); it!=local.end(); ++it)
            values[it->first]=it->second;
        loading.erase(canonical);                  // 本文件加载完毕，允许别的分支再引用它（只禁止环）
    }

    std::map<std::string,std::string> values;       ///< 合并后的参数表（不含环境变量覆盖）
    std::map<std::string,bool> allowedKeys;         ///< defaults 登记的合法参数名
};

} // namespace water_fangsha

#endif // WATER_FANGSHA_PARAMETER_CONFIG_H
