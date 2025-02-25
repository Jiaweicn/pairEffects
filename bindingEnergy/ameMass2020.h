//created by DeepSeek online service
//processing AME2020 data

// 核素数据结构
struct Nuclide2020 {
   int A;
   int Z;
   string element;
   double mass_excess;
   double binding_energy_per_nucleon;
};

// 解析单行数据
bool parseLine2020(const string& line, Nuclide2020& result) {
   // 固定列格式解析
   try {
      // 提取关键字段（基于文件格式描述）
      //int NZ = stoi(line.substr(2, 2));// 列2-4: NZ
      //int N = stoi(line.substr(6, 3));// 列5-9: N
      int Z = stoi(line.substr(11, 3));// 列10-14: Z
      int A = stoi(line.substr(16, 3));// 列15-19: A
      string el = line.substr(20, 2);// 列21-23: 元素符号
      string binding_str = line.substr(55, 13);// binding energy per nucleon(keV)

      // 清理元素符号中的空格
      el.erase(remove(el.begin(), el.end(), ' '), el.end());
      if (el.empty()) return false;

      // 处理Binding energy中的特殊字符(如#和*)
      size_t invalid_pos = binding_str.find_first_of("#*");
      if (invalid_pos != string::npos)
          binding_str = binding_str.substr(0, invalid_pos);
      if (binding_str.empty()) return false;
      double bindingEnergyPerNucleon = stod(binding_str);// 转换为数值

      // 填充结果
      result.A = A;
      result.Z = Z;
      result.element = el;
      result.binding_energy_per_nucleon = bindingEnergyPerNucleon;
      return true;
   } catch (...) {
      return false; // 解析失败
   }
}

// 主函数
bool getNucleusBindingE2020(short a, short z,string &eleName,double &bindingEnergyPerNucleon){
   const char massFile[60]="/home/cai/prjs/pairEffect/ameMass2020.txt";//https://www-nds.iaea.org/amdc/ame2020/mass_1.mas20.txt
   ifstream file(massFile);
   if (!file) {
      cerr << "Error: Cannot open "<<massFile <<"\n";
      return false;
   }

   string line;
   bool in_data_section = false;

   // 跳过文件头，定位到数据区域
   while (getline(file, line)) {
      if (line.find("MASS EXCESS") != string::npos) {
         in_data_section = true;
         getline(file, line); // 跳过列标题行
         break;
      }
   }

   if (!in_data_section) return false;

   // 逐行搜索目标核素
   while (getline(file, line)){
      Nuclide2020 nuclide;
      if (parseLine2020(line, nuclide)) {
         if (nuclide.A == a && nuclide.Z == z) {
             eleName=nuclide.element;
             bindingEnergyPerNucleon=nuclide.binding_energy_per_nucleon*1e-3;//MeV
             return true;
           }
       }
   }
   return false;
}
