import fmt
import k8s.io.apimachinery.pkg.runtime.schema as schema

class KindVisitor:
    def visit_daemon_set(self, kind):
        pass
    
    def visit_deployment(self, kind):
        pass
    
    def visit_job(self, kind):
        pass
    
    def visit_pod(self, kind):
        pass
    
    def visit_replica_set(self, kind):
        pass
    
    def visit_replication_controller(self, kind):
        pass
    
    def visit_stateful_set(self, kind):
        pass
    
    def visit_cron_job(self, kind):
        pass

class GroupKindElement:
    def __init__(self, group, kind):
        self.group = group
        self.kind = kind

    def accept(self, visitor):
        if (self.group_match("apps", "extensions") and self.kind == "DaemonSet"):
            visitor.visit_daemon_set(self)
        elif (self.group_match("apps", "extensions") and self.kind == "Deployment"):
            visitor.visit_deployment(self)
        elif (self.group_match("batch") and self.kind == "Job"):
            visitor.visit_job(self)
        elif (self.group_match("", "core") and self.kind == "Pod"):
            visitor.visit_pod(self)
        elif (self.group_match("apps", "extensions") and self.kind == "ReplicaSet"):
            visitor.visit_replica_set(self)
        elif (self.group_match("", "core") and self.kind == "ReplicationController"):
            visitor.visit_replication_controller(self)
        elif (self.group_match("apps") and self.kind == "StatefulSet"):
            visitor.visit_stateful_set(self)
        elif (self.group_match("batch") and self.kind == "CronJob"):
            visitor.visit_cron_job(self)
        else:
            return ValueError("no visitor method exists for {}".format(self))
        return None
    
    def group_match(self, *groups):
        for g in groups:
            if self.group == g:
                return True
        return False
