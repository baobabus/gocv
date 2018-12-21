#include "features2d.h"

AKAZE AKAZE_Create() {
    // TODO: params
    return new cv::Ptr<cv::AKAZE>(cv::AKAZE::create());
}

void AKAZE_Close(AKAZE a) {
    delete a;
}

struct KeyPoints AKAZE_Detect(AKAZE a, Mat src) {
    std::vector<cv::KeyPoint> detected;
    (*a)->detect(*src, detected);

    KeyPoint* kps = new KeyPoint[detected.size()];

    for (size_t i = 0; i < detected.size(); ++i) {
        KeyPoint k = {detected[i].pt.x, detected[i].pt.y, detected[i].size, detected[i].angle,
                      detected[i].response, detected[i].octave, detected[i].class_id
                     };
        kps[i] = k;
    }

    KeyPoints ret = {kps, (int)detected.size()};
    return ret;
}

struct KeyPoints AKAZE_DetectAndCompute(AKAZE a, Mat src, Mat mask, Mat desc) {
    std::vector<cv::KeyPoint> detected;
    (*a)->detectAndCompute(*src, *mask, detected, *desc);

    KeyPoint* kps = new KeyPoint[detected.size()];

    for (size_t i = 0; i < detected.size(); ++i) {
        KeyPoint k = {detected[i].pt.x, detected[i].pt.y, detected[i].size, detected[i].angle,
                      detected[i].response, detected[i].octave, detected[i].class_id
                     };
        kps[i] = k;
    }

    KeyPoints ret = {kps, (int)detected.size()};
    return ret;
}

AgastFeatureDetector AgastFeatureDetector_Create() {
    // TODO: params
    return new cv::Ptr<cv::AgastFeatureDetector>(cv::AgastFeatureDetector::create());
}

void AgastFeatureDetector_Close(AgastFeatureDetector a) {
    delete a;
}

struct KeyPoints AgastFeatureDetector_Detect(AgastFeatureDetector a, Mat src) {
    std::vector<cv::KeyPoint> detected;
    (*a)->detect(*src, detected);

    KeyPoint* kps = new KeyPoint[detected.size()];

    for (size_t i = 0; i < detected.size(); ++i) {
        KeyPoint k = {detected[i].pt.x, detected[i].pt.y, detected[i].size, detected[i].angle,
                      detected[i].response, detected[i].octave, detected[i].class_id
                     };
        kps[i] = k;
    }

    KeyPoints ret = {kps, (int)detected.size()};
    return ret;
}

BRISK BRISK_Create() {
    // TODO: params
    return new cv::Ptr<cv::BRISK>(cv::BRISK::create());
}

void BRISK_Close(BRISK b) {
    delete b;
}

struct KeyPoints BRISK_Detect(BRISK b, Mat src) {
    std::vector<cv::KeyPoint> detected;
    (*b)->detect(*src, detected);

    KeyPoint* kps = new KeyPoint[detected.size()];

    for (size_t i = 0; i < detected.size(); ++i) {
        KeyPoint k = {detected[i].pt.x, detected[i].pt.y, detected[i].size, detected[i].angle,
                      detected[i].response, detected[i].octave, detected[i].class_id
                     };
        kps[i] = k;
    }

    KeyPoints ret = {kps, (int)detected.size()};
    return ret;
}

struct KeyPoints BRISK_DetectAndCompute(BRISK b, Mat src, Mat mask, Mat desc) {
    std::vector<cv::KeyPoint> detected;
    (*b)->detectAndCompute(*src, *mask, detected, *desc);

    KeyPoint* kps = new KeyPoint[detected.size()];

    for (size_t i = 0; i < detected.size(); ++i) {
        KeyPoint k = {detected[i].pt.x, detected[i].pt.y, detected[i].size, detected[i].angle,
                      detected[i].response, detected[i].octave, detected[i].class_id
                     };
        kps[i] = k;
    }

    KeyPoints ret = {kps, (int)detected.size()};
    return ret;
}

GFTTDetector GFTTDetector_Create() {
    // TODO: params
    return new cv::Ptr<cv::GFTTDetector>(cv::GFTTDetector::create());
}

void GFTTDetector_Close(GFTTDetector a) {
    delete a;
}

struct KeyPoints GFTTDetector_Detect(GFTTDetector a, Mat src) {
    std::vector<cv::KeyPoint> detected;
    (*a)->detect(*src, detected);

    KeyPoint* kps = new KeyPoint[detected.size()];

    for (size_t i = 0; i < detected.size(); ++i) {
        KeyPoint k = {detected[i].pt.x, detected[i].pt.y, detected[i].size, detected[i].angle,
                      detected[i].response, detected[i].octave, detected[i].class_id
                     };
        kps[i] = k;
    }

    KeyPoints ret = {kps, (int)detected.size()};
    return ret;
}

KAZE KAZE_Create() {
    // TODO: params
    return new cv::Ptr<cv::KAZE>(cv::KAZE::create());
}

void KAZE_Close(KAZE a) {
    delete a;
}

struct KeyPoints KAZE_Detect(KAZE a, Mat src) {
    std::vector<cv::KeyPoint> detected;
    (*a)->detect(*src, detected);

    KeyPoint* kps = new KeyPoint[detected.size()];

    for (size_t i = 0; i < detected.size(); ++i) {
        KeyPoint k = {detected[i].pt.x, detected[i].pt.y, detected[i].size, detected[i].angle,
                      detected[i].response, detected[i].octave, detected[i].class_id
                     };
        kps[i] = k;
    }

    KeyPoints ret = {kps, (int)detected.size()};
    return ret;
}

struct KeyPoints KAZE_DetectAndCompute(KAZE a, Mat src, Mat mask, Mat desc) {
    std::vector<cv::KeyPoint> detected;
    (*a)->detectAndCompute(*src, *mask, detected, *desc);

    KeyPoint* kps = new KeyPoint[detected.size()];

    for (size_t i = 0; i < detected.size(); ++i) {
        KeyPoint k = {detected[i].pt.x, detected[i].pt.y, detected[i].size, detected[i].angle,
                      detected[i].response, detected[i].octave, detected[i].class_id
                     };
        kps[i] = k;
    }

    KeyPoints ret = {kps, (int)detected.size()};
    return ret;
}

MSER MSER_Create() {
    // TODO: params
    return new cv::Ptr<cv::MSER>(cv::MSER::create());
}

void MSER_Close(MSER a) {
    delete a;
}

struct KeyPoints MSER_Detect(MSER a, Mat src) {
    std::vector<cv::KeyPoint> detected;
    (*a)->detect(*src, detected);

    KeyPoint* kps = new KeyPoint[detected.size()];

    for (size_t i = 0; i < detected.size(); ++i) {
        KeyPoint k = {detected[i].pt.x, detected[i].pt.y, detected[i].size, detected[i].angle,
                      detected[i].response, detected[i].octave, detected[i].class_id
                     };
        kps[i] = k;
    }

    KeyPoints ret = {kps, (int)detected.size()};
    return ret;
}

FastFeatureDetector FastFeatureDetector_Create() {
    // TODO: params
    return new cv::Ptr<cv::FastFeatureDetector>(cv::FastFeatureDetector::create());
}

void FastFeatureDetector_Close(FastFeatureDetector f) {
    delete f;
}

struct KeyPoints FastFeatureDetector_Detect(FastFeatureDetector f, Mat src) {
    std::vector<cv::KeyPoint> detected;
    (*f)->detect(*src, detected);

    KeyPoint* kps = new KeyPoint[detected.size()];

    for (size_t i = 0; i < detected.size(); ++i) {
        KeyPoint k = {detected[i].pt.x, detected[i].pt.y, detected[i].size, detected[i].angle,
                      detected[i].response, detected[i].octave, detected[i].class_id
                     };
        kps[i] = k;
    }

    KeyPoints ret = {kps, (int)detected.size()};
    return ret;
}

ORB ORB_Create() {
    // TODO: params
    return new cv::Ptr<cv::ORB>(cv::ORB::create());
}

void ORB_Close(ORB o) {
    delete o;
}

struct KeyPoints ORB_Detect(ORB o, Mat src) {
    std::vector<cv::KeyPoint> detected;
    (*o)->detect(*src, detected);

    KeyPoint* kps = new KeyPoint[detected.size()];

    for (size_t i = 0; i < detected.size(); ++i) {
        KeyPoint k = {detected[i].pt.x, detected[i].pt.y, detected[i].size, detected[i].angle,
                      detected[i].response, detected[i].octave, detected[i].class_id
                     };
        kps[i] = k;
    }

    KeyPoints ret = {kps, (int)detected.size()};
    return ret;
}

struct KeyPoints ORB_DetectAndCompute(ORB o, Mat src, Mat mask, Mat desc) {
    std::vector<cv::KeyPoint> detected;
    (*o)->detectAndCompute(*src, *mask, detected, *desc);

    KeyPoint* kps = new KeyPoint[detected.size()];

    for (size_t i = 0; i < detected.size(); ++i) {
        KeyPoint k = {detected[i].pt.x, detected[i].pt.y, detected[i].size, detected[i].angle,
                      detected[i].response, detected[i].octave, detected[i].class_id
                     };
        kps[i] = k;
    }

    KeyPoints ret = {kps, (int)detected.size()};
    return ret;
}

SimpleBlobDetectorParams *_defaultSimpleBlobDetectorParams;

SimpleBlobDetectorParams *defaultSimpleBlobDetectorParams() {
    // not reentrant!
    if (_defaultSimpleBlobDetectorParams == NULL) {
        cv::SimpleBlobDetector::Params params;
        _defaultSimpleBlobDetectorParams = new SimpleBlobDetectorParams();
        _defaultSimpleBlobDetectorParams->thresholdStep = params.thresholdStep;
        _defaultSimpleBlobDetectorParams->minThreshold = params.minThreshold;
        _defaultSimpleBlobDetectorParams->maxThreshold = params.maxThreshold;
        _defaultSimpleBlobDetectorParams->minRepeatability = params.minRepeatability;
        _defaultSimpleBlobDetectorParams->minDistBetweenBlobs = params.minDistBetweenBlobs;
        _defaultSimpleBlobDetectorParams->filterByColor = params.filterByColor;
        _defaultSimpleBlobDetectorParams->blobColor = params.blobColor;
        _defaultSimpleBlobDetectorParams->filterByArea = params.filterByArea;
        _defaultSimpleBlobDetectorParams->minArea = params.minArea;
        _defaultSimpleBlobDetectorParams->maxArea = params.maxArea;
        _defaultSimpleBlobDetectorParams->filterByCircularity = params.filterByCircularity;
        _defaultSimpleBlobDetectorParams->minCircularity = params.minCircularity;
        _defaultSimpleBlobDetectorParams->maxCircularity = params.maxCircularity;
        _defaultSimpleBlobDetectorParams->filterByInertia = params.filterByInertia;
        _defaultSimpleBlobDetectorParams->minInertiaRatio = params.minInertiaRatio;
        _defaultSimpleBlobDetectorParams->maxInertiaRatio = params.maxInertiaRatio;
        _defaultSimpleBlobDetectorParams->filterByConvexity = params.filterByConvexity;
        _defaultSimpleBlobDetectorParams->minConvexity = params.minConvexity;
        _defaultSimpleBlobDetectorParams->maxConvexity = params.maxConvexity;
    }
    return _defaultSimpleBlobDetectorParams;
}

SimpleBlobDetector SimpleBlobDetector_Create(SimpleBlobDetectorParams *params) {
    // TODO: params
    cv::SimpleBlobDetector::Params pars;
    if (params) {
        pars.thresholdStep = params->thresholdStep;
        pars.minThreshold = params->minThreshold;
        pars.maxThreshold = params->maxThreshold;
        pars.minRepeatability = params->minRepeatability;
        pars.minDistBetweenBlobs = params->minDistBetweenBlobs;
        pars.filterByColor = params->filterByColor;
        pars.blobColor = params->blobColor;
        pars.filterByArea = params->filterByArea;
        pars.minArea = params->minArea;
        pars.maxArea = params->maxArea;
        pars.filterByCircularity = params->filterByCircularity;
        pars.minCircularity = params->minCircularity;
        pars.maxCircularity = params->maxCircularity;
        pars.filterByInertia = params->filterByInertia;
        pars.minInertiaRatio = params->minInertiaRatio;
        pars.maxInertiaRatio = params->maxInertiaRatio;
        pars.filterByConvexity = params->filterByConvexity;
        pars.minConvexity = params->minConvexity;
        pars.maxConvexity = params->maxConvexity;
    }
    return new cv::Ptr<cv::SimpleBlobDetector>(cv::SimpleBlobDetector::create(pars));
}

void SimpleBlobDetector_Close(SimpleBlobDetector b) {
    delete b;
}

struct KeyPoints SimpleBlobDetector_Detect(SimpleBlobDetector b, Mat src) {
    std::vector<cv::KeyPoint> detected;
    (*b)->detect(*src, detected);

    KeyPoint* kps = new KeyPoint[detected.size()];

    for (size_t i = 0; i < detected.size(); ++i) {
        KeyPoint k = {detected[i].pt.x, detected[i].pt.y, detected[i].size, detected[i].angle,
                      detected[i].response, detected[i].octave, detected[i].class_id
                     };
        kps[i] = k;
    }

    KeyPoints ret = {kps, (int)detected.size()};
    return ret;
}

BFMatcher BFMatcher_Create() {
    return new cv::Ptr<cv::BFMatcher>(cv::BFMatcher::create());
}

BFMatcher BFMatcher_CreateWithParams(int normType, bool crossCheck) {
    return new cv::Ptr<cv::BFMatcher>(cv::BFMatcher::create(normType, crossCheck));
}

void BFMatcher_Close(BFMatcher b) {
    delete b;
}

struct MultiDMatches BFMatcher_KnnMatch(BFMatcher b, Mat query, Mat train, int k) {
    std::vector< std::vector<cv::DMatch> > matches;
    (*b)->knnMatch(*query, *train, matches, k);

    DMatches *dms = new DMatches[matches.size()];
    for (size_t i = 0; i < matches.size(); ++i) {
        DMatch *dmatches = new DMatch[matches[i].size()];
        for (size_t j = 0; j < matches[i].size(); ++j) {
            DMatch dmatch = {matches[i][j].queryIdx, matches[i][j].trainIdx, matches[i][j].imgIdx,
                             matches[i][j].distance};
            dmatches[j] = dmatch;
        }
        dms[i] = {dmatches, (int) matches[i].size()};
    }
    MultiDMatches ret = {dms, (int) matches.size()};
    return ret;
}

struct MultiDMatches BFMatcher_KnnMatchWithParams(BFMatcher b, Mat query, Mat train, int k, Mat mask, bool compactResult) {
    std::vector< std::vector<cv::DMatch> > matches;
    (*b)->knnMatch(*query, *train, matches, k, *mask, compactResult);

    DMatches *dms = new DMatches[matches.size()];
    for (size_t i = 0; i < matches.size(); ++i) {
        DMatch *dmatches = new DMatch[matches[i].size()];
        for (size_t j = 0; j < matches[i].size(); ++j) {
            DMatch dmatch = {matches[i][j].queryIdx, matches[i][j].trainIdx, matches[i][j].imgIdx,
                             matches[i][j].distance};
            dmatches[j] = dmatch;
        }
        dms[i] = {dmatches, (int) matches[i].size()};
    }
    MultiDMatches ret = {dms, (int) matches.size()};
    return ret;
}

void DrawKeyPoints(Mat src, struct KeyPoints kp, Mat dst, Scalar s, int flags) {
        std::vector<cv::KeyPoint> keypts;
        cv::KeyPoint keypt;

        for (int i = 0; i < kp.length; ++i) {
                keypt = cv::KeyPoint(kp.keypoints[i].x, kp.keypoints[i].y,
                                kp.keypoints[i].size, kp.keypoints[i].angle, kp.keypoints[i].response,
                                kp.keypoints[i].octave, kp.keypoints[i].classID);
                keypts.push_back(keypt);
        }

        cv::Scalar color = cv::Scalar(s.val1, s.val2, s.val3, s.val4);

        cv::drawKeypoints(*src, keypts, *dst, color, static_cast<cv::DrawMatchesFlags>(flags));
}
