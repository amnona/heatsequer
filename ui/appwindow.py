<?xml version="1.0" encoding="UTF-8"?>
<ui version="4.0">
 <class>MainWindow</class>
 <widget class="QMainWindow" name="MainWindow">
  <property name="geometry">
   <rect>
    <x>0</x>
    <y>0</y>
    <width>800</width>
    <height>600</height>
   </rect>
  </property>
  <property name="windowTitle">
   <string>MainWindow</string>
  </property>
  <widget class="QWidget" name="centralwidget">
   <widget class="QLabel" name="label">
    <property name="geometry">
     <rect>
      <x>30</x>
      <y>40</y>
      <width>91</width>
      <height>16</height>
     </rect>
    </property>
    <property name="text">
     <string>Experiments:</string>
    </property>
   </widget>
   <widget class="QPushButton" name="bMainLoadNew">
    <property name="geometry">
     <rect>
      <x>170</x>
      <y>30</y>
      <width>115</width>
      <height>32</height>
     </rect>
    </property>
    <property name="text">
     <string>Load New</string>
    </property>
   </widget>
   <widget class="QPushButton" name="bMainPlot">
    <property name="geometry">
     <rect>
      <x>30</x>
      <y>260</y>
      <width>131</width>
      <height>51</height>
     </rect>
    </property>
    <property name="text">
     <string>Plot</string>
    </property>
   </widget>
   <widget class="QPushButton" name="bMainAdvancedPlot">
    <property name="geometry">
     <rect>
      <x>160</x>
      <y>260</y>
      <width>131</width>
      <height>51</height>
     </rect>
    </property>
    <property name="text">
     <string>Advanced Plot</string>
    </property>
   </widget>
   <widget class="QListWidget" name="bMainList">
    <property name="geometry">
     <rect>
      <x>30</x>
      <y>60</y>
      <width>256</width>
      <height>192</height>
     </rect>
    </property>
    <property name="selectionMode">
     <enum>QAbstractItemView::ExtendedSelection</enum>
    </property>
   </widget>
   <widget class="QPushButton" name="bPickleLoad">
    <property name="geometry">
     <rect>
      <x>170</x>
      <y>0</y>
      <width>115</width>
      <height>32</height>
     </rect>
    </property>
    <property name="text">
     <string>PickleLoad</string>
    </property>
   </widget>
   <widget class="QTabWidget" name="tabWidget">
    <property name="geometry">
     <rect>
      <x>330</x>
      <y>40</y>
      <width>181</width>
      <height>251</height>
     </rect>
    </property>
    <property name="currentIndex">
     <number>0</number>
    </property>
    <widget class="QWidget" name="tab">
     <attribute name="title">
      <string>Samples</string>
     </attribute>
     <widget class="QPushButton" name="bMainFilterSamples">
      <property name="geometry">
       <rect>
        <x>0</x>
        <y>0</y>
        <width>121</width>
        <height>32</height>
       </rect>
      </property>
      <property name="text">
       <string>Filter Samples</string>
      </property>
     </widget>
     <widget class="QPushButton" name="bMainSortSamples">
      <property name="geometry">
       <rect>
        <x>0</x>
        <y>30</y>
        <width>121</width>
        <height>32</height>
       </rect>
      </property>
      <property name="text">
       <string>Sort Samples</string>
      </property>
     </widget>
     <widget class="QPushButton" name="bJoinFields">
      <property name="geometry">
       <rect>
        <x>0</x>
        <y>60</y>
        <width>115</width>
        <height>32</height>
       </rect>
      </property>
      <property name="text">
       <string>Join Fields</string>
      </property>
     </widget>
     <widget class="QPushButton" name="bFilterOrigReads">
      <property name="geometry">
       <rect>
        <x>0</x>
        <y>90</y>
        <width>121</width>
        <height>32</height>
       </rect>
      </property>
      <property name="text">
       <string>Filter OrigReads</string>
      </property>
     </widget>
    </widget>
    <widget class="QWidget" name="tab_2">
     <attribute name="title">
      <string>Bacteria</string>
     </attribute>
     <widget class="QPushButton" name="bMainClusterBacteria">
      <property name="geometry">
       <rect>
        <x>10</x>
        <y>10</y>
        <width>121</width>
        <height>32</height>
       </rect>
      </property>
      <property name="text">
       <string>Cluster Bacteria</string>
      </property>
     </widget>
     <widget class="QPushButton" name="bMainFilterMinReads">
      <property name="geometry">
       <rect>
        <x>10</x>
        <y>40</y>
        <width>121</width>
        <height>32</height>
       </rect>
      </property>
      <property name="text">
       <string>Filter MinReads</string>
      </property>
     </widget>
     <widget class="QPushButton" name="bMainFilterTaxonomy">
      <property name="geometry">
       <rect>
        <x>10</x>
        <y>70</y>
        <width>121</width>
        <height>32</height>
       </rect>
      </property>
      <property name="text">
       <string>Filter Taxonomy</string>
      </property>
     </widget>
     <widget class="QPushButton" name="bFilterPresence">
      <property name="geometry">
       <rect>
        <x>10</x>
        <y>100</y>
        <width>121</width>
        <height>32</height>
       </rect>
      </property>
      <property name="text">
       <string>Filter Presence</string>
      </property>
     </widget>
     <widget class="QPushButton" name="bFilterMean">
      <property name="geometry">
       <rect>
        <x>10</x>
        <y>130</y>
        <width>121</width>
        <height>32</height>
       </rect>
      </property>
      <property name="text">
       <string>Filter Mean</string>
      </property>
     </widget>
    </widget>
   </widget>
  </widget>
  <widget class="QMenuBar" name="menubar">
   <property name="geometry">
    <rect>
     <x>0</x>
     <y>0</y>
     <width>800</width>
     <height>22</height>
    </rect>
   </property>
  </widget>
  <widget class="QStatusBar" name="statusbar"/>
 </widget>
 <resources/>
 <connections/>
</ui>
